# run_texas.jl — full study on ACTIVSg2000: SQ under three queue orders, CO (i)-(iii)
#   MODELS=sq,co1,co2,co3  (default all)   RUN_TAG=name of the results folder
#   quick test:  N_LOADS=20 TOTAL_GW=8 HOUR_STRIDE=6 julia --project=. run_texas.jl

include("batch0.jl")

const MODELS  = split(get(ENV, "MODELS", "sq,co1,co2,co3"), ',')
const RUN_TAG = get(ENV, "RUN_TAG", "texas_$(LOADS_SOURCE)_n$(N_LOADS)_s$(HOUR_STRIDE)_seed$(SEED)")

net   = load_network()
sites = LOADS_SOURCE == "dc" ? load_dc_sites(net) : nothing
sc    = build_scenario(net; sites=sites)
lds   = LOADS_SOURCE == "dc" ? dc_loads(net, sites) : sample_loads(net)
dir = mkpath(joinpath(OUT_DIR, RUN_TAG))
CSV.write(joinpath(dir, "loads.csv"), DataFrame([(k=l.k, bus=net.bus_id[l.bus], area=l.area, Lmax=l.Lmax,
          L2028=l.L[1], L2030=l.L[2], L2032=l.L[3], pclr=l.pclr, phi=l.phi, theta=l.theta) for l in lds]))

rows = NamedTuple[]
Yall = 1:length(sc.years)

# ---- status quo, three queue orders ----
if "sq" in MODELS
    for rule in (:random, :smallest, :largest)
        t0 = time()
        order = queue_order(lds, rule)
        a_sq, a_raw, amax = run_sq(net, sc, lds, order; tag="SQ-$rule")
        sol = evaluate_fixed(net, sc, lds, a_sq; tag="SQ-$rule eval")
        strand = stranded(net, sc, lds, a_sq, a_raw, amax, order; tag="SQ-$rule stranded")
        append!(rows, summary_rows("SQ_$rule", sol, sc, lds; amax=amax, stranded=strand, solve_s=time() - t0))
        write_allocations(dir, "SQ_$rule", sol, sc, lds; amax=amax)
        CSV.write(joinpath(dir, "summary.csv"), DataFrame(rows))          # rewrite after each model, keeps partial results
    end
end

# ---- coordinated model, nested variants ----
prev, prev_var = nothing, :none
brows = NamedTuple[]
# WARM_FROM=<results dir>/alloc_CO_i.csv warm-starts the first CO variant of this job from a variant
# solved in an earlier job, so CO (i), (ii), (iii) can run as separate jobs that each fit in mit_normal
if haskey(ENV, "WARM_FROM")
    prev = read_alloc_start(ENV["WARM_FROM"], sc, lds)
    prev_var = Symbol(basename(ENV["WARM_FROM"]))
end
for (key, var) in (("co1", :i), ("co2", :ii), ("co3", :iii))
    key in MODELS || continue
    t0 = time()
    m = build_co(net, sc, lds; variant=var)
    # warm start from the previous CO variant in this job: its allocation is feasible here (f = 0 and
    # x = 0 remain allowed), so Gurobi always has a solution to improve and a time limit never leaves nothing
    if prev !== nothing
        for k in eachindex(lds), yi in Yall
            set_start_value(m[:a][k, yi], prev.a[k, yi]); set_start_value(m[:f][k, yi], prev.f[k, yi])
            set_start_value(m[:z][k, yi], prev.z[k, yi])
        end
        @printf("[CO-%s] warm start from CO-%s\n", var, prev_var)
    end
    st = solve_n1!(m, net, sc, Yall; tag="CO-$var", relax_first=true)
    rv = get(m.obj_dict, :report_vals, nothing)          # set when the result is the start's saved values
    if !has_values(m) && rv === nothing       # e.g. TIME_LIMIT before any solution: record it and go on
        @printf("[CO-%s] status %s with no solution after %.0fs; row left NaN\n", var, st, time() - t0)
        nanv = fill(NaN, length(Yall))
        push!(rows, [(model="CO_$var", year=y, firm_GW=NaN, total_GW=NaN, curt_served_avg_GW=NaN, MW_equiv_GW=NaN,
                      n_admitted=0, n_zeroed_by_threshold=0, stranded_GW=NaN, cost_M=NaN, welfare_M=NaN,
                      served_value_M=NaN, gen_cost_M=NaN, upg_cost_M=NaN, upg_GW=NaN, shed_MWh=NaN,
                      solve_s=time() - t0) for y in sc.years]...)
        CSV.write(joinpath(dir, "summary.csv"), DataFrame(rows))
        continue
    end
    val = rv === nothing ? value : (v -> rv[v])
    sol = read_solution(m, net, sc, lds; val=val)
    obj = rv === nothing ? objective_value(m) : value(val, objective_function(m))
    @printf("[CO-%s] status %s, obj %.4e, %.0fs\n", var, st, obj, time() - t0)
    # optimality certificate: LP relaxation bound, value of the constructed start, final value
    ub = get(m.obj_dict, :lp_bound, NaN)
    # Gurobi's best bound after branching is the one to report; the LP bound ignores the threshold and is weak
    bb = get(m.obj_dict, :mip_bound, NaN)
    isnan(bb) && (bb = try objective_bound(m) catch; NaN end)
    get(m.obj_dict, :used_start, false) && @printf("[CO-%s] result is the constructed start (MILP did not improve it)\n", var)
    push!(brows, (model="CO_$var", status=string(st), lp_bound=ub, start_obj=get(m.obj_dict, :start_obj, NaN),
                  obj=obj, mip_bound=bb, mip_gap_pct=100 * (bb - obj) / max(abs(obj), 1.0),
                  gap_to_lp_pct=100 * (ub - obj) / max(abs(ub), 1.0)))
    CSV.write(joinpath(dir, "bounds.csv"), DataFrame(brows))
    append!(rows, summary_rows("CO_$var", sol, sc, lds; solve_s=time() - t0))
    duals = (var == :iii && rv === nothing) ? cap_duals(m, net, sc, lds, sol.z; tag="CO-iii duals") : nothing
    write_allocations(dir, "CO_$var", sol, sc, lds; duals=duals)
    CSV.write(joinpath(dir, "summary.csv"), DataFrame(rows))
    global prev, prev_var = (a=sol.a, f=sol.f, z=sol.z), var
    m = nothing; GC.gc()                      # free this model before building the next one
end

println("\n", DataFrame(rows))
println("results in $dir")
