# metrics.jl — welfare evaluation, summary tables, CSV output

# Pull the solution of a solved CO model into plain arrays. val reads a variable's value: JuMP's value for
# a solved model, or a lookup into saved values (the constructed start, if the model has no solution).
function read_solution(m, net::Network, sc::Scenario, lds; val=value)
    ev(x::VariableRef) = val(x)
    ev(x::Real) = Float64(x)
    ev(x) = value(val, x)                        # affine expression
    K, Y = 1:length(lds), 1:length(sc.years)
    a = [ev(m[:a][k, yi]) for k in K, yi in Y]
    f = [ev(m[:f][k, yi]) for k in K, yi in Y]
    z = [round(ev(m[:z][k, yi])) for k in K, yi in Y]
    qavg = [sum(sc.Nt[t] / N_YR * ev(m[:q][k, yi, t]) for t in 1:sc.T) for k in K, yi in Y]   # avg curtailable MW served
    x = [ev(m[:x][l, yi]) for l in 1:net.nl, yi in Y]
    yrs = Dict(:served => [ev(m[:served][yi]) for yi in Y],
               :gencost => [ev(m[:gencost][yi]) for yi in Y],
               :upgcost => [ev(m[:upgcost][yi]) for yi in Y],
               :shed_mwh => [sum(sc.Nt[t] * ev(m[:shed][n, yi, t]) for n in 1:net.nb, t in 1:sc.T) for yi in Y])
    return (a=a, f=f, z=z, qavg=qavg, x=x, yrs=yrs)
end

# Welfare of the SQ allocation: CO (i) with a and z fixed, so only dispatch is chosen (an LP).
# If lowering an earlier year broke feasibility, it shows up as shed.
function evaluate_fixed(net::Network, sc::Scenario, lds, a_fix; tag="eval")
    m = build_co(net, sc, lds; variant=:i, shed=true)     # SQ's allocation may need a little shed; keep the slack
    fix_allocation!(m; z=Float64.(a_fix .> 0), a=a_fix)
    st = solve_n1!(m, net, sc, 1:length(sc.years); tag=tag)
    st == MOI.OPTIMAL || @printf("[%s] status %s\n", tag, st)
    has_values(m) && return read_solution(m, net, sc, lds)
    @printf("[%s] WARNING: no solution, this model's row will be NaN\n", tag)   # keep the other models running
    K, nY = length(lds), length(sc.years)
    nanv = fill(NaN, nY)
    return (a=a_fix, f=zeros(K, nY), z=Float64.(a_fix .> 0), qavg=zeros(K, nY), x=zeros(net.nl, nY),
            yrs=Dict(:served => nanv, :gencost => nanv, :upgcost => nanv, :shed_mwh => nanv))
end

# Marginal value of headroom: fix z at the CO optimum, solve the LP, read the dual of (6).
# (11) links the years, so a single year's dual is not unique; compare loads on the sum over years.
function cap_duals(m, net::Network, sc::Scenario, lds, z_opt; tag="duals")
    fix_allocation!(m; z=z_opt)
    solve_n1!(m, net, sc, 1:length(sc.years); tag=tag)
    return [shadow_price(m[:cap][k, yi]) for k in 1:length(lds), yi in 1:length(sc.years)]   # $/MW of request
end

# Stranded MW: after SQ, the most MW the loads zeroed by the threshold could still take
# together, with every SQ allocation fixed and N-1 enforced. One LP per year. Summing each zeroed
# load's own offer would double count headroom the zeroed loads share.
stranded(net, sc, lds, a_sq, a_raw, amax, order; tag="stranded") =
    SQ_RULE == "lp" ? stranded_mw(net, sc, lds, a_sq, amax; tag=tag) :
                      stranded_prop(net, sc, lds, a_raw, amax, order; tag=tag)

function stranded_mw(net::Network, sc::Scenario, lds, a_sq, amax; tag="stranded")
    K, nY = length(lds), length(sc.years)
    out = zeros(nY)
    for yi in 1:nY
        Z = [k for k in 1:K if !lds[k].pclr && 1e-3 < amax[k, yi] < lds[k].theta]
        isempty(Z) && continue
        m = new_model(; method=SQ_METHOD)
        @variable(m, a[k=1:K] >= 0)
        for k in 1:K
            k in Z ? set_upper_bound(a[k], lds[k].L[yi]) : fix(a[k], a_sq[k, yi]; force=true)
        end
        loads_at = [Int[] for _ in 1:net.nb]; for k in 1:K; push!(loads_at[lds[k].bus], k); end
        wd(n, y, t) = isempty(loads_at[n]) ? AffExpr(0.0) : sum(a[k] for k in loads_at[n])
        add_network!(m, net, sc, yi:yi, wd)
        for l in 1:net.nl; fix(m[:x][l, yi], 0.0; force=true); end
        shed_avg = @expression(m, sum(sc.Nt[t] / N_YR * m[:shed][n, yi, t] for n in 1:net.nb, t in 1:sc.T))
        @objective(m, Max, sum(a[k] for k in Z) - SQ_SHED_PEN * shed_avg)
        st = solve_n1!(m, net, sc, yi:yi; tag="$tag y$(sc.years[yi])")
        out[yi] = st == MOI.OPTIMAL ? sum(value(a[k]) for k in Z) : NaN
    end
    return out
end

# One row per (model, year)
# Welfare is linear in v_k when v_k is equal for all loads, so MW_equiv and cost_M are reported too:
#   welfare(v) = v * MW_equiv - cost_M   when v_k is equal for all loads,
# and any v can be checked without rerunning. MW_equiv = firm MW + beta * average curtailable MW served.
function summary_rows(label, sol, sc::Scenario, lds; amax=nothing, stranded=nothing, solve_s=NaN)
    rows = NamedTuple[]
    for (yi, y) in enumerate(sc.years)
        zeroed = amax === nothing ? 0 : sum(!lds[k].pclr && 1e-3 < amax[k, yi] < lds[k].theta for k in eachindex(lds))
        w = sol.yrs[:served][yi] - sol.yrs[:gencost][yi] - sol.yrs[:upgcost][yi]
        push!(rows, (model=label, year=y,
                     firm_GW=sum(sol.a[:, yi]) / 1e3,
                     total_GW=sum(sol.a[:, yi] .+ sol.f[:, yi]) / 1e3,
                     curt_served_avg_GW=sum(sol.qavg[:, yi]) / 1e3,
                     MW_equiv_GW=(sum(sol.a[:, yi]) + BETA * sum(sol.qavg[:, yi])) / 1e3,
                     n_admitted=sum(sol.a[:, yi] .+ sol.f[:, yi] .> 1e-3),
                     n_zeroed_by_threshold=zeroed,
                     stranded_GW=stranded === nothing ? 0.0 : stranded[yi] / 1e3,
                     cost_M=(sol.yrs[:gencost][yi] + sol.yrs[:upgcost][yi]) / 1e6,
                     welfare_M=w / 1e6,
                     served_value_M=sol.yrs[:served][yi] / 1e6,
                     gen_cost_M=sol.yrs[:gencost][yi] / 1e6,
                     upg_cost_M=sol.yrs[:upgcost][yi] / 1e6,
                     upg_GW=sum(sol.x[:, yi]) / 1e3,
                     shed_MWh=sol.yrs[:shed_mwh][yi],
                     solve_s=solve_s))
    end
    return rows
end

# Per-load allocations. With equal v_k, CO's optimum is usually not unique in WHICH loads get the MW
# so these files are one optimal solution; compare models on the totals in summary.csv.
function write_allocations(dir, label, sol, sc::Scenario, lds; amax=nothing, duals=nothing)
    rows = NamedTuple[]
    for ld in lds, (yi, y) in enumerate(sc.years)
        k = ld.k
        push!(rows, (k=k, bus=ld.bus, area=ld.area, pclr=ld.pclr, year=y, L=ld.L[yi], theta=ld.theta,
                     a=sol.a[k, yi], f=sol.f[k, yi], q_avg=sol.qavg[k, yi],
                     amax=amax === nothing ? NaN : amax[k, yi], cap_dual=duals === nothing ? NaN : duals[k, yi]))
    end
    CSV.write(joinpath(dir, "alloc_$(label).csv"), DataFrame(rows))
end

# Read a previous run's alloc_<model>.csv back as a warm start (a, f, z), for WARM_FROM in run_texas.jl.
# Same data, seed and stride give the same loads, so k and year line up.
function read_alloc_start(path, sc::Scenario, lds)
    df = CSV.read(path, DataFrame)
    K, nY = length(lds), length(sc.years)
    a = zeros(K, nY); f = zeros(K, nY)
    for r in eachrow(df)
        yi = findfirst(==(r.year), sc.years)
        yi === nothing && continue
        a[r.k, yi] = r.a; f[r.k, yi] = r.f
    end
    @assert maximum(df.k) == K "WARM_FROM has $(maximum(df.k)) loads, this run has $K"
    return (a=a, f=f, z=Float64.(a .+ f .> 1e-3))
end
