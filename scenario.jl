# scenario.jl — representative hours, background demand, gen availability, large-load requests

using CSV, DataFrames, Random

struct Scenario
    years::Vector{Int}
    T::Int
    Nt::Vector{Float64}              # N_t, real hours each modeled hour stands for, sums to 8760
    stress::Vector{Bool}             # t \in \mathcal S
    D::Array{Float64,3}              # D[b, yi, t] background demand [MW]
    Pbar::Matrix{Float64}            # \bar P[g, t] [MW] (same every study year, no new gen modeled)
    day::Vector{String}; hour::Vector{Int}
end

struct LargeLoad
    k::Int
    bus::Int
    area::Int
    L::Vector{Float64}               # L_{k,y} per study year [MW]
    Lmax::Float64                    # L_k
    pclr::Bool                       # k \in \mathcal K^P
    theta::Float64                   # \theta_k
    v::Float64                       # v_k
end

theta_k(Lmax) = Lmax > THRESH_MW ? THRESH_MW : THRESH_SHR * Lmax

# ACTIVSg2000 areas are ERCOT's weather zones
const AREA_ZONE = Dict(1 => "FWES", 2 => "NRTH", 3 => "WEST", 4 => "SOUT",
                       5 => "NCEN", 6 => "SCEN", 7 => "COAS", 8 => "EAST")

# Data center sites (dc_sites.csv) placed on the case and scaled to it.
# The synthetic case's zones are smaller than ERCOT's real 2025 zones (e.g. Far West 1.3 GW in the
# case vs 8.8 GW in 2025), so every site's MW is multiplied by its zone's factor
#   f_z = (case load in area z) / (2025 peak of zone z)
# which keeps each site the same share of its zone. Sites go to the nearest bus of at least
# MIN_KV_LOAD in their zone if COORD_CSV exists, otherwise to a random such bus in the zone.
function load_dc_sites(net::Network; seed=SEED)
    st = CSV.read(DC_CSV, DataFrame)
    zp = Dict(r.zone => r.peak_mw for r in eachrow(CSV.read(ZPEAK_CSV, DataFrame)))
    zone_area = Dict(v => k for (k, v) in AREA_ZONE)
    area_pd = Dict(a => sum(max.(net.Pd[net.area .== a], 0.0)) for a in keys(AREA_ZONE))
    fz = Dict(z => area_pd[zone_area[z]] / zp[z] for z in keys(zone_area))

    cand = Dict(a => findall(i -> net.area[i] == a && net.kv[i] >= MIN_KV_LOAD, 1:net.nb) for a in keys(AREA_ZONE))
    coords = isfile(COORD_CSV) ? CSV.read(COORD_CSV, DataFrame) : nothing
    if coords !== nothing
        idx = Dict(b => i for (i, b) in enumerate(net.bus_id))
        lat = fill(NaN, net.nb); lon = fill(NaN, net.nb)
        for r in eachrow(coords); haskey(idx, r.bus_id) && (lat[idx[r.bus_id]] = r.lat; lon[idx[r.bus_id]] = r.lon); end
    end
    hav(la1, lo1, la2, lo2) = 2 * 6371 * asin(sqrt(sind((la2 - la1) / 2)^2 + cosd(la1) * cosd(la2) * sind((lo2 - lo1) / 2)^2))

    rng = MersenneTwister(seed + 7)
    st.area = [zone_area[z] for z in st.zone]
    st.scale = [fz[z] for z in st.zone]
    st.bus = map(eachrow(st)) do r
        c = cand[r.area]
        coords === nothing && return rand(rng, c)
        c[argmin([isnan(lat[i]) ? Inf : hav(r.lat, r.lon, lat[i], lon[i]) for i in c])]
    end
    st.existing_case = st.existing_mw .* st.scale
    st.planned_case = st.planned_mw .* st.scale
    @printf("[DC ] %d sites placed %s; zone scale factors %s\n", nrow(st),
            coords === nothing ? "at RANDOM buses in their zone (no COORD_CSV)" : "at the nearest bus",
            join(["$z=$(round(fz[z], digits=2))" for z in sort(collect(keys(fz)))], " "))
    @printf("[DC ] existing %.1f GW real -> %.1f GW in the case; planned %.1f GW real -> %.1f GW\n",
            sum(st.existing_mw) / 1e3, sum(st.existing_case) / 1e3, sum(st.planned_mw) / 1e3, sum(st.planned_case) / 1e3)
    return st
end

function build_scenario(net::Network; sites=nothing, profiles=PROFILE_CSV, years=YEARS, stride=HOUR_STRIDE)
    pr = CSV.read(profiles, DataFrame)       # all representative hours; thinned below, after stress hours are known
    T  = nrow(pr)
    Nt = Float64.(pr.days_weight)            # one modeled hour per real hour of each represented day
    @assert abs(sum(Nt) - N_YR) < 1e-6 "N_t must sum to 8760, got $(sum(Nt))"

    # background demand. Each bus follows its own zone's 2025 hourly shape (load / zone peak), with the
    # case's Pd taken as the zone peak, grown from BASE_YEAR. Existing data centers are base loads:
    # their (scaled) MW is placed at their buses, flat in every hour, and the same share of zone load
    # is taken out of the other buses, so each zone's total is unchanged at its peak.
    zonal = all(c -> c in names(pr), ["lf_" * z for z in values(AREA_ZONE)])
    lf(b, t) = zonal && haskey(AREA_ZONE, net.area[b]) ? pr[t, "lf_" * AREA_ZONE[net.area[b]]] : pr.load_factor[t]
    exist_bus = zeros(net.nb); share = Dict(a => 0.0 for a in keys(AREA_ZONE))
    if sites !== nothing && EXISTING_AT_SITES
        for r in eachrow(sites); exist_bus[r.bus] += r.existing_case; end
        for a in keys(AREA_ZONE)
            pd = sum(max.(net.Pd[net.area .== a], 0.0))
            share[a] = pd > 0 ? sum(exist_bus[net.area .== a]) / pd : 0.0
        end
    end
    D = zeros(net.nb, length(years), T)
    for (yi, y) in enumerate(years), t in 1:T, b in 1:net.nb
        g = (1 + GROWTH)^(y - BASE_YEAR)
        D[b, yi, t] = max(net.Pd[b], 0.0) * (1 - get(share, net.area[b], 0.0)) * lf(b, t) * g + exist_bus[b]
    end

    # available output: wind/solar follow the profile, hydro derated, thermal at Pmax
    Pbar = zeros(net.ng, T)
    for g in 1:net.ng, t in 1:T
        f = net.fuel[g]
        Pbar[g, t] = net.Pmax[g] * (f == "wind" ? pr.wind_cf[t] : f == "solar" ? pr.solar_cf[t] :
                                    f == "hydro" ? HYDRO_AF : 1.0)
    end

    # stress hours: top STRESS_FRAC of net load (demand minus wind and solar), first study year
    netload = [sum(D[:, 1, t]) - sum(Pbar[g, t] for g in 1:net.ng if net.fuel[g] in ("wind", "solar"); init=0.0) for t in 1:T]
    nS = max(1, round(Int, STRESS_FRAC * T))
    stress = falses(T); stress[sortperm(netload, rev=true)[1:nS]] .= true

    # thin hours (HOUR_STRIDE > 1): keep every stride-th hour of each day AND every stress hour, and give
    # each dropped hour's weight to the nearest kept hour of the same day, so N_t still sums to 8760 and
    # the peak hours that bind N-1 are never dropped
    if stride > 1
        keep = [((pr.hour[t] - 1) % stride == 0) || stress[t] for t in 1:T]
        w = zeros(T)
        for t in 1:T
            same = [u for u in 1:T if keep[u] && pr.day[u] == pr.day[t]]
            u = same[argmin([abs(pr.hour[v] - pr.hour[t]) for v in same])]
            w[u] += Nt[t]
        end
        idx = findall(keep)
        pr, D, Pbar, stress, Nt = pr[idx, :], D[:, :, idx], Pbar[:, idx], stress[idx], w[idx]
        T = length(idx)
        @assert abs(sum(Nt) - N_YR) < 1e-6
    end

    @printf("[SCN] %d hours/yr (stride %d), %d stress hours, peak background demand %.1f GW in %d%s\n",
            T, stride, nS, maximum(sum(D[:, end, :], dims=1)) / 1e3, years[end],
            zonal ? ", zonal 2025 shapes" : ", one system-wide shape")
    sites === nothing || @printf("[SCN] existing data centers: %s\n", EXISTING_AT_SITES ?
        @sprintf("%.1f GW placed at their buses as base loads", sum(exist_bus) / 1e3) : "left inside zone load (EXISTING_AT_SITES=false)")
    flush(stdout)
    return Scenario(years, T, Nt, stress, D, Pbar, String.(pr.day), pr.hour)
end

# Sample N requests: zone by share, bus uniform among high-voltage buses in the zone,
# size lognormal then scaled so the 2032 total is TOTAL_GW, ramp nondecreasing over study years
function sample_loads(net::Network; n=N_LOADS, total_gw=TOTAL_GW, seed=SEED, zone_csv=ZONE_CSV)
    rng = MersenneTwister(seed)
    zs = CSV.read(zone_csv, DataFrame)
    cand = Dict(a => findall(i -> net.area[i] == a && net.kv[i] >= MIN_KV_LOAD, 1:net.nb) for a in zs.area)
    w = cumsum(zs.share ./ sum(zs.share))

    sz = clamp.(exp.(log(300.0) .+ 0.8 .* randn(rng, n)), 75.0, 1500.0)    # MW, ERCOT large load >= 75 MW
    sz .*= total_gw * 1e3 / sum(sz)

    loads = LargeLoad[]
    for k in 1:n
        a = zs.area[findfirst(>=(rand(rng)), w)]
        bus = rand(rng, cand[a])
        r1 = 0.2 + 0.4 * rand(rng); r2 = r1 + (1 - r1) * rand(rng)          # share of request by 2028, 2030
        L = sz[k] .* [r1, r2, 1.0]
        pclr = rand(rng) < P_PCLR
        push!(loads, LargeLoad(k, bus, a, L, sz[k], pclr, theta_k(sz[k]), V_FIRM))
    end
    @printf("[LDS] %d requests, %.1f GW in %d (%.1f GW in %d), %d PCLR\n", n, sum(l.L[end] for l in loads) / 1e3,
            YEARS[end], sum(l.L[1] for l in loads) / 1e3, YEARS[1], sum(l.pclr for l in loads))
    flush(stdout)
    return loads
end

# Studied loads from the data center sites: every site with at least MIN_SITE_MW planned (real MW),
# at its scaled size. No in-service dates are available, so the 2028 and 2030 shares of the request
# are drawn at random (nondecreasing), as is the PCLR election. N_LOADS > 0 keeps a random subset.
function dc_loads(net::Network, sites; n=N_LOADS, seed=SEED)
    rng = MersenneTwister(seed)
    pick = findall(sites.planned_mw .>= MIN_SITE_MW)
    n > 0 && n < length(pick) && (pick = sort(shuffle(rng, pick)[1:n]))
    loads = LargeLoad[]
    for (k, i) in enumerate(pick)
        r = sites[i, :]
        sz = r.planned_case
        r1 = 0.2 + 0.4 * rand(rng); r2 = r1 + (1 - r1) * rand(rng)
        pclr = rand(rng) < P_PCLR
        push!(loads, LargeLoad(k, r.bus, r.area, sz .* [r1, r2, 1.0], sz, pclr, theta_k(sz), V_FIRM))
    end
    @printf("[LDS] %d data center requests, %.1f GW in %d (%.1f GW in %d), %d PCLR, %d below the %.0f MW threshold size\n",
            length(loads), sum(l.L[end] for l in loads) / 1e3, YEARS[end], sum(l.L[1] for l in loads) / 1e3, YEARS[1],
            sum(l.pclr for l in loads), sum(l.Lmax <= THRESH_MW for l in loads), THRESH_MW)
    flush(stdout)
    return loads
end
