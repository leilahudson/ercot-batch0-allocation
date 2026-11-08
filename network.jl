# network.jl — MATPOWER reader, DC shift factors (PTDF) and outage factors (LODF)

using LinearAlgebra, SparseArrays, Statistics

struct Network
    nb::Int; nl::Int; ng::Int
    bus_id::Vector{Int}          # MATPOWER bus number
    area::Vector{Int}            # ACTIVSg2000 area = ERCOT weather zone (1 Far West ... 8 East)
    kv::Vector{Float64}
    Pd::Vector{Float64}          # MW, case snapshot
    slack::Int
    fr::Vector{Int}; to::Vector{Int}     # line endpoints (bus index)
    b::Vector{Float64}                   # 1/x [MW/rad on 100 MVA base]
    rate::Vector{Float64}                # \bar F_l [MW]
    expand::Vector{Bool}                 # l \in \mathcal L^exp
    C::Vector{Float64}                   # C_l [$/MW-yr]
    xbar::Vector{Float64}                # \bar x_l [MW]
    gbus::Vector{Int}; Pmax::Vector{Float64}; gamma::Vector{Float64}; fuel::Vector{String}
    PTDF::Matrix{Float64}                # nl x nb
    LODF::Matrix{Float64}                # nl x nl, column c = outage of line c
    cont_ok::Vector{Bool}                # false if outage of c islands a bus (radial), skipped in N-1
end

# Pull one numeric block, e.g. mpc.bus = [ ... ];
function _mp_block(txt, name)
    m = match(Regex("mpc\\.$name\\s*=\\s*\\[(.*?)\\];", "s"), txt)
    rows = Vector{Vector{Float64}}()
    for ln in split(m.captures[1], '\n')
        s = strip(rstrip(strip(split(ln, '%')[1]), ';'))
        isempty(s) || push!(rows, parse.(Float64, split(s)))
    end
    return permutedims(reduce(hcat, rows))      # rows x cols
end

function _mp_cells(txt, name)
    m = match(Regex("mpc\\.$name\\s*=\\s*\\{(.*?)\\};", "s"), txt)
    m === nothing && return String[]
    return [strip(strip(s), [';', '\'']) for s in split(m.captures[1]) if !isempty(strip(s))]
end

function load_network(file::String=CASE_FILE)
    txt = read(file, String)
    base = parse(Float64, match(r"mpc\.baseMVA\s*=\s*([\d.]+)", txt).captures[1])
    B = _mp_block(txt, "bus"); G = _mp_block(txt, "gen"); L = _mp_block(txt, "branch")
    gc = _mp_block(txt, "gencost"); fuel = String.(_mp_cells(txt, "genfuel"))
    isempty(fuel) && (fuel = fill("unk", size(G, 1)))

    # buses
    bus_id = Int.(B[:, 1]); idx = Dict(b => i for (i, b) in enumerate(bus_id))
    slack = findfirst(==(3), Int.(B[:, 2]))

    # in-service branches only
    on = L[:, 11] .> 0; L = L[on, :]
    fr = [idx[Int(v)] for v in L[:, 1]]; to = [idx[Int(v)] for v in L[:, 2]]
    tap = [t == 0 ? 1.0 : t for t in L[:, 9]]
    b = base ./ (L[:, 4] .* tap)                 # MW per rad, phase shifts ignored
    rate = L[:, 6]
    expand = L[:, 9] .== 0                       # lines can be upgraded, transformers cannot
    xmed = median(L[expand, 4])
    C = UPG_COST .* L[:, 4] ./ xmed              # reactance as a length proxy (no line lengths in the case)
    xbar = UPG_MAX .* rate .* expand

    # generators: all units, including ones off in the snapshot (planning view)
    gbus = [idx[Int(v)] for v in G[:, 1]]
    Pmax = G[:, 9]
    # gencost model 2, ncost 3: [c2 c1 c0]; gamma = average marginal cost over [0, Pmax]
    gamma = [gc[i, 1] == 2 ? gc[i, 6] + gc[i, 5] * Pmax[i] : 0.0 for i in 1:size(G, 1)]

    nb, nl = length(bus_id), length(fr)
    PTDF, LODF, cont_ok = shift_factors(nb, nl, fr, to, b, slack)

    @printf("[NET] %d buses, %d branches (%d expandable), %d gens, %.1f GW load, %.1f GW gen, %d radial outages skipped\n",
            nb, nl, sum(expand), length(gbus), sum(B[:, 3]) / 1e3, sum(Pmax) / 1e3, sum(.!cont_ok))
    flush(stdout)
    return Network(nb, nl, length(gbus), bus_id, Int.(B[:, 7]), B[:, 10], B[:, 3], slack,
                   fr, to, b, rate, expand, C, xbar, gbus, Pmax, gamma, fuel, PTDF, LODF, cont_ok)
end

# PTDF = Bf * Bbus^{-1} with the slack column zero; LODF_lc = H_lc / (1 - H_cc), H = PTDF * A'
function shift_factors(nb, nl, fr, to, b, slack)
    A  = sparse(1:nl, fr, 1.0, nl, nb) - sparse(1:nl, to, 1.0, nl, nb)   # incidence, +1 at from bus
    Bf = spdiagm(b) * A
    Bbus = A' * Bf
    ns = setdiff(1:nb, slack)
    PTDF = zeros(nl, nb)
    PTDF[:, ns] = (lu(Bbus[ns, ns]) \ Matrix(Bf[:, ns]'))'          # Bbus symmetric
    H = PTDF * A'                                # H[l,c] = PTDF[l,fr_c] - PTDF[l,to_c]
    den = 1 .- diag(H)
    cont_ok = den .> 1e-5                        # den ~ 0 means line c is a bridge (radial)
    LODF = H ./ max.(den, 1e-5)'
    LODF[:, .!cont_ok] .= 0.0
    LODF[diagind(LODF)] .= 0.0                   # l = c never checked
    return PTDF, LODF, cont_ok
end

# Small network built in code (toy example, Sec. 5)
function make_network(; kv, area, Pd, fr, to, x, rate, gbus, Pmax, gamma, fuel=fill("ng", length(gbus)),
                        slack=1, C=fill(UPG_COST, length(fr)), xbar=UPG_MAX .* rate)
    nb, nl = length(Pd), length(fr)
    b = 100.0 ./ x
    PTDF, LODF, cont_ok = shift_factors(nb, nl, fr, to, b, slack)
    return Network(nb, nl, length(gbus), collect(1:nb), area, kv, Pd, slack, fr, to, b, rate,
                   fill(true, nl), C, xbar, gbus, Pmax, gamma, fuel, PTDF, LODF, cont_ok)
end
