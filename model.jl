# model.jl — DC network core, N-1 constraint generation, and the coordinated allocation model CO

using JuMP, HiGHS, Printf

function new_model(; method::Int=GRB_METHOD)
    if SOLVER == "gurobi"
        m = Model(() -> Gurobi.Optimizer(GRB_ENV[]))
        set_optimizer_attribute(m, "Threads", THREADS)
        set_optimizer_attribute(m, "MIPGap", MIP_GAP)
        set_optimizer_attribute(m, "TimeLimit", TIME_LIMIT)
        set_optimizer_attribute(m, "Method", method)
        GRB_NODEFILE > 0 && set_optimizer_attribute(m, "NodefileStart", GRB_NODEFILE)   # caps B&B tree memory
    elseif SOLVER == "glpk"
        m = Model(GLPK.Optimizer)
        set_optimizer_attribute(m, "mip_gap", MIP_GAP)
        set_optimizer_attribute(m, "tm_lim", round(Int, 1000 * TIME_LIMIT))
    else
        m = Model(HiGHS.Optimizer)
        set_optimizer_attribute(m, "threads", THREADS)
        set_optimizer_attribute(m, "mip_rel_gap", MIP_GAP)
        set_optimizer_attribute(m, "time_limit", TIME_LIMIT)
    end
    VERBOSE || set_silent(m)
    return m
end

# Network core for study years Y (indices into sc.years).
# wd(n, yi, t) returns the large-load withdrawal at bus n as an AffExpr (0 if none).
# Flows use the angle form F_l = b_l(θ_fr - θ_to), which equals sum_b PTDF_lb π_b
# once the nodal balances hold, so (2) and the PTDF flow definition are both enforced.
function add_network!(m, net::Network, sc::Scenario, Y, wd; shed::Bool=true)
    T = sc.T
    @variable(m, θ[1:net.nb, Y, 1:T])
    @variable(m, 0 <= p[g=1:net.ng, yi=Y, t=1:T] <= sc.Pbar[g, t])                   # (3)
    if shed      # VOLL slack, not in the proposal; 2000 x years x hours variables, so CO leaves it out
        @variable(m, 0 <= shed[n=1:net.nb, yi=Y, t=1:T] <= sc.D[n, yi, t])
    else
        m[:shed] = JuMP.Containers.DenseAxisArray([AffExpr(0.0) for n in 1:net.nb, yi in Y, t in 1:T], 1:net.nb, Y, 1:T)
    end
    shed = m[:shed]
    @variable(m, 0 <= x[l=1:net.nl, yi=Y] <= net.xbar[l])                            # (12) upper bound
    for yi in Y, t in 1:T
        fix(θ[net.slack, yi, t], 0.0; force=true)
    end

    @expression(m, F[l=1:net.nl, yi=Y, t=1:T], net.b[l] * (θ[net.fr[l], yi, t] - θ[net.to[l], yi, t]))

    # nodal balance: generation - demand + shed - large load = net flow out
    gens_at = [Int[] for _ in 1:net.nb]; for g in 1:net.ng; push!(gens_at[net.gbus[g]], g); end
    out_l   = [Int[] for _ in 1:net.nb]; in_l = [Int[] for _ in 1:net.nb]
    for l in 1:net.nl; push!(out_l[net.fr[l]], l); push!(in_l[net.to[l]], l); end
    @constraint(m, bal[n=1:net.nb, yi=Y, t=1:T],
        sum(p[g, yi, t] for g in gens_at[n]; init=0.0) - sc.D[n, yi, t] + shed[n, yi, t] - wd(n, yi, t)
        == sum(F[l, yi, t] for l in out_l[n]; init=0.0) - sum(F[l, yi, t] for l in in_l[n]; init=0.0))

    # (4) intact limits, both directions
    @constraint(m, flo_up[l=1:net.nl, yi=Y, t=1:T],  F[l, yi, t] <= net.rate[l] + x[l, yi])
    @constraint(m, flo_dn[l=1:net.nl, yi=Y, t=1:T], -F[l, yi, t] <= net.rate[l] + x[l, yi])

    # (12) upgrades stay once built; non-expandable lines get none
    Yv = collect(Y)
    for i in 1:length(Yv)-1, l in 1:net.nl
        @constraint(m, x[l, Yv[i]] <= x[l, Yv[i+1]])
    end
    for l in 1:net.nl
        net.expand[l] || for yi in Y; fix(x[l, yi], 0.0; force=true); end
    end

    m[:n1_keys] = Set{NTuple{5,Int}}()          # (l, c, yi, t, sign) already added
    return m
end

# One round of N-1 checking: post-contingency flow F_l + LODF_lc F_c against α(\bar F_l + x_l), eq. (5).
# For each (line, year, hour) only the worst outage and only the violated side get a cut.
# val(v) reads a variable's value (value, or callback_value inside a lazy callback); add(con) adds the
# cut (a model constraint, or a lazy constraint). Returns (cuts added, largest violation in MW).
function find_n1_cuts!(m, net::Network, sc::Scenario, Y, val, add; keys=m[:n1_keys])
    F, x, θ = m[:F], m[:x], m[:θ]
    cand = Tuple{Float64,Int,Int,Int,Int,Int}[]       # (violation, l, c, yi, t, sign)
    Fv = zeros(net.nl)
    for yi in Y
        lim = ALPHA .* (net.rate .+ [val(x[l, yi]) for l in 1:net.nl])
        for t in 1:sc.T
            th = [val(θ[n, yi, t]) for n in 1:net.nb]
            @. Fv = net.b * (th[net.fr] - th[net.to])     # same as the value of F[l, yi, t], much faster
            worst_v = fill(N1_TOL, net.nl); worst_c = zeros(Int, net.nl); worst_s = zeros(Int, net.nl)
            for c in 1:net.nl
                (net.cont_ok[c] && abs(Fv[c]) > 1e-6) || continue
                fc = Fv[c]
                @inbounds for l in 1:net.nl
                    post = Fv[l] + net.LODF[l, c] * fc
                    v = abs(post) - lim[l]
                    if v > worst_v[l]
                        worst_v[l] = v; worst_c[l] = c; worst_s[l] = post > 0 ? 1 : -1
                    end
                end
            end
            for l in 1:net.nl
                worst_c[l] > 0 && push!(cand, (worst_v[l], l, worst_c[l], yi, t, worst_s[l]))
            end
        end
    end
    sort!(cand, rev=true)
    added = 0; vmax = isempty(cand) ? 0.0 : cand[1][1]
    for (v, l, c, yi, t, s) in cand
        added >= N1_MAXCUTS && break
        key = (l, c, yi, t, s)
        key in keys && continue                       # already in the model, violation is solver tolerance
        # cut carries the same N1_TOL the check accepts, so a solution accepted earlier stays feasible
        # when this cut is added later (matters in SQ, where earlier loads are fixed at their LP values)
        add(@build_constraint(s * (F[l, yi, t] + net.LODF[l, c] * F[c, yi, t]) <= ALPHA * (net.rate[l] + x[l, yi]) + N1_TOL))
        push!(keys, key); added += 1
    end
    return added, vmax
end

add_n1_cuts!(m, net::Network, sc::Scenario, Y) = find_n1_cuts!(m, net, sc, Y, value, con -> add_constraint(m, con))

# MILP stage with lazy N-1 cuts: one branch and bound; each integer solution Gurobi finds is checked
# and, if it violates N-1, cut off with the violated limits. The final incumbent is therefore N-1 secure,
# also when the solve stops at TIME_LIMIT. Cuts from the LP rounds stay in the model as ordinary rows.
function solve_milp_lazy!(m, net::Network, sc::Scenario, Y; tag="")
    active = Ref(true); ncb = Ref(0); ncut = Ref(0)
    lazy_keys = Set{NTuple{5,Int}}()  # lazy cuts are not rows of the JuMP model, so keep them out of m[:n1_keys]
    function cb(cb_data)
        active[] || return
        callback_node_status(cb_data, m) == MOI.CALLBACK_NODE_STATUS_INTEGER || return
        ncb[] += 1
        added, vmax = find_n1_cuts!(m, net, sc, Y, v -> callback_value(cb_data, v),
                                    con -> MOI.submit(m, MOI.LazyConstraint(cb_data), con); keys=lazy_keys)
        ncut[] += added
        N1_LOG && added > 0 && @printf("[N-1 %s] lazy: incumbent %d cut off by %d cuts (max viol %.1f MW)\n", tag, ncb[], added, vmax)
        flush(stdout)
    end
    MOI.set(m, MOI.LazyConstraintCallback(), cb)
    set_time_limit_sec(m, MILP_TIME)
    t0 = time(); optimize!(m)
    active[] = false                  # later solves of this model (duals) are LPs; keep the callback silent
    st = termination_status(m)
    @printf("[N-1 %s] MILP with lazy N-1: status %s, %d integer solutions checked, %d lazy cuts, %.0fs\n",
            tag, st, ncb[], ncut[], time() - t0)
    has_values(m) && @printf("[N-1 %s] incumbent %.4e, best bound %.4e, gap %.2f%%\n", tag, objective_value(m),
                             objective_bound(m), 100 * relative_gap(m))
    if has_values(m)                  # check the final solution against every outage once more (adds nothing)
        added, vmax = find_n1_cuts!(m, net, sc, Y, value, con -> nothing; keys=Set{NTuple{5,Int}}())
        vmax > 0.1 && @printf("[N-1 %s] WARNING: final solution violates N-1 by %.1f MW (%d limits)\n", tag, vmax, added)   # below 0.1 MW is solver tolerance
    end
    return st
end

# Solve, check N-1, add cuts, repeat. With relax_first, the cut rounds run on the LP relaxation
# first (cheap), then integrality is restored and rounds continue on the MILP.
function solve_n1!(m, net::Network, sc::Scenario, Y; tag="", relax_first=false)
    set_time_limit_sec(m, TIME_LIMIT)                  # the MILP stage may have shortened it
    undo = relax_first ? relax_integrality(m) : nothing
    for r in 1:N1_ROUNDS
        flush(stdout)                                  # SLURM logs are block-buffered
        t0 = time(); optimize!(m)
        st = termination_status(m)
        if !(st in (MOI.OPTIMAL, MOI.TIME_LIMIT) && has_values(m))
            @printf("[N-1 %s] round %d: status %s, stopping\n", tag, r, st); return st
        end
        obj = objective_value(m)                       # read before cuts modify the model
        t_solve = time() - t0; t_grb = solve_time(m); t1 = time()
        added, vmax = N1_ON ? add_n1_cuts!(m, net, sc, Y) : (0, 0.0)
        N1_LOG && @printf("[N-1 %s] round %d: obj %.4e, %d cuts added (max viol %.1f MW), solve %.1fs (solver %.1fs), check %.1fs\n",
                           tag, r, obj, added, vmax, t_solve, t_grb, time() - t1)
        if added == 0
            if undo !== nothing
                # LP rounds done: obj is the relaxation's bound. Read the LP allocation before the model changes,
                # then build a feasible incumbent for branch and bound
                lpa = haskey(m, :z) ? (value.(m[:a]), value.(m[:f])) : nothing
                undo(); undo = nothing                # back to MILP, keep the cuts
                lpa === nothing || round_and_repair!(m, net, sc, Y, obj, lpa...; tag=tag)
                if N1_LAZY && SOLVER in ("gurobi", "glpk")
                    st = solve_milp_lazy!(m, net, sc, Y; tag=tag)
                    return use_start_if_better!(m, net, sc, Y, st; tag=tag)
                end
                continue
            end
            return st
        end
    end
    undo === nothing || undo()
    @printf("[N-1 %s] hit N1_ROUNDS=%d with violations left\n", tag, N1_ROUNDS)
    optimize!(m)
    return termination_status(m)
end

# N-1 cut rounds on the model as it is (used with z fixed). Returns (status, converged), where converged
# means the last solve had values and the N-1 check found nothing to add: the solution is then feasible
# even if that solve stopped at the time limit.
function cut_rounds!(m, net::Network, sc::Scenario, Y; tag="")
    st = MOI.OPTIMIZE_NOT_CALLED
    for r in 1:N1_ROUNDS
        flush(stdout); t0 = time(); optimize!(m)
        st = termination_status(m)
        (st in (MOI.OPTIMAL, MOI.TIME_LIMIT) && has_values(m)) || return st, false
        obj = objective_value(m)
        added, vmax = add_n1_cuts!(m, net, sc, Y)
        N1_LOG && @printf("[N-1 %s] round %d: obj %.4e, %d cuts added (max viol %.1f MW), %.1fs\n", tag, r, obj, added, vmax, time() - t0)
        added == 0 && return st, true
    end
    return MOI.OTHER_LIMIT, false
end

# Feasible start for CO's MILP. Candidates for z, tried in order: (1) the LP relaxation rounded, z = 1 where a load already gets at
# least its threshold (any MW for a PCLR); (2) the warm start from the previous variant, if set;
# (3) z = 0, which admits no load and is feasible whenever background demand is N-1 secure.
# Each candidate is solved with z fixed (an LP) and N-1 cut rounds; the best one becomes the MIP start.
# Its cuts stay in the model, so the start satisfies every N-1 limit found so far.
# The LP bound and this start's value give a proven optimality gap even if the MILP does not improve.
function round_and_repair!(m, net::Network, sc::Scenario, Y, lp_bound, aLP, fLP; tag="")
    lds = m[:lds]; z = m[:z]
    K, nY = length(lds), length(Y)
    aL = [aLP[k, yi] for k in 1:K, yi in Y]; fL = [fLP[k, yi] for k in 1:K, yi in Y]
    function rounded(frac)            # z = 1 where the LP gives at least frac * threshold (any MW for a PCLR)
        zr = [lds[k].pclr ? Float64(aL[k, yi] + fL[k, yi] > 1e-3) : Float64(aL[k, yi] >= frac * lds[k].theta - 1e-6)
              for k in 1:K, yi in 1:nY]
        for k in 1:K, yi in 2:nY; zr[k, yi] = max(zr[k, yi], zr[k, yi-1]); end       # (11)
        return zr
    end
    # the relaxation spreads MW thinly (it ignores the threshold), so also try admitting loads at half the
    # threshold; the z-fixed LP then decides whether they fit
    # order: the previous variant's solution first (feasible here, since f = 0 and x = 0 stay allowed), then the
    # LP roundings, then no loads as a last resort
    cands = Pair{String,Matrix{Float64}}[]
    zs = start_value(z[1, first(Y)])
    zs === nothing || push!(cands, "warm start" => [round(something(start_value(z[k, yi]), 0.0)) for k in 1:K, yi in Y])
    push!(cands, "LP rounded" => rounded(1.0), "LP rounded at half threshold" => rounded(0.5), "no loads" => zeros(K, nY))

    vars = all_variables(m)
    set_time_limit_sec(m, START_TIME)                 # each candidate solve; a failed candidate is skipped
    best_obj, best_vals, best_name = -Inf, nothing, ""
    for (name, zc) in cands
        best_vals !== nothing && name == "no loads" && break       # only needed if nothing else worked
        for k in 1:K, (j, yi) in enumerate(Y); fix(z[k, yi], zc[k, j]; force=true); end
        # with z fixed the model is an LP; relax the (fixed) binaries so Gurobi solves it as one, not as a MIP
        undo_c = relax_integrality(m)
        st, conv = cut_rounds!(m, net, sc, Y; tag="$tag start ($name)")
        ok = conv && has_values(m)                     # feasible: N-1 check passed on a solution with values
        obj = ok ? objective_value(m) : -Inf
        # keep only the decisions; Gurobi completes dispatch and angles itself. A start copying all ~10^5
        # continuous values was rejected over a 1e-5 MW balance residual (above its 1e-6 tolerance).
        ok && (vals_c = Dict(s => value.(m[s]) for s in (:z, :a, :f, :x)))
        ok && (full_c = value.(vars))                  # every variable, to report this start if all else fails
        undo_c()
        @printf("[N-1 %s] start candidate '%s': %d of %d load-years admitted, status %s, obj %.4e\n",
                tag, name, round(Int, sum(zc)), K * nY, st, obj)
        if ok && obj > best_obj
            best_obj, best_name = obj, name
            best_vals = vals_c
            m[:start_z] = copy(zc)
            m[:start_vals] = Dict(zip(vars, full_c))
        end
    end
    for k in 1:K, yi in Y; unfix(z[k, yi]); end                     # binary again, no bounds needed
    m[:lp_bound] = lp_bound
    if best_vals === nothing
        @printf("[N-1 %s] WARNING: no feasible start found\n", tag)
        m[:start_obj] = NaN
    else
        if MOI.supports(backend(m), MOI.VariablePrimalStart(), MOI.VariableIndex)
            for s in (:z, :a, :f, :x), i in eachindex(m[s])
                v = best_vals[s][i]
                set_start_value(m[s][i], s == :z ? round(v) : max(v, 0.0))
            end
        end
        m[:start_obj] = best_obj
        @printf("[N-1 %s] MIP start from '%s', obj %.4e; LP bound %.4e, gap %.2f%%\n", tag, best_name, best_obj,
                lp_bound, 100 * (lp_bound - best_obj) / max(abs(lp_bound), 1.0))
    end
    flush(stdout)
end

# If branch and bound returned nothing, or something worse than the constructed start (Gurobi may decline
# a MIP start), re-solve with the start's z fixed (an LP; its N-1 cuts are already rows of the model) and
# report that solution. Its value and the bounds in bounds.csv still certify how far it is from optimal.
function use_start_if_better!(m, net::Network, sc::Scenario, Y, st; tag="")
    m[:mip_bound] = try objective_bound(m) catch; NaN end          # Gurobi's bound from the MILP attempt
    zs = get(m.obj_dict, :start_z, nothing)
    so = get(m.obj_dict, :start_obj, NaN)
    zs === nothing && return st
    if has_values(m) && objective_value(m) >= so - 1e-7 * abs(so)
        return st                                                    # the MILP did at least as well
    end
    @printf("[N-1 %s] MILP returned %s; using the constructed start (obj %.4e)\n", tag,
            has_values(m) ? @sprintf("a worse solution (%.4e)", objective_value(m)) : "no solution", so)
    z = m[:z]
    for k in axes(z, 1), (j, yi) in enumerate(Y); fix(z[k, yi], zs[k, j]; force=true); end
    relax_integrality(m)             # z is fixed, so this is the start's LP; left relaxed so the values stay valid
    set_time_limit_sec(m, TIME_LIMIT)
    set_time_limit_sec(m, START_TIME)   # the start was already solved once; don't let a slow re-solve eat the job
    st2, conv = cut_rounds!(m, net, sc, Y; tag="$tag start re-solve")
    m[:used_start] = true
    if !(conv && has_values(m))
        # report the start from its saved values; they passed every N-1 check when it was built
        @printf("[N-1 %s] start re-solve status %s; reporting the start from its saved values\n", tag, st2)
        m[:report_vals] = m[:start_vals]
    end
    return MOI.TIME_LIMIT                                            # feasible, not proven optimal
end

# Coordinated model, eq. (1)-(12).
#   variant :i   firm only   (f = 0, x = 0)
#   variant :ii  + curtailable (x = 0)
#   variant :iii + upgrades
function build_co(net::Network, sc::Scenario, lds::Vector{LargeLoad}; variant::Symbol=:iii, shed::Bool=CO_SHED)
    m = new_model()
    K, Y, T = 1:length(lds), 1:length(sc.years), 1:sc.T
    @variable(m, a[K, Y] >= 0)                       # firm allocation
    @variable(m, f[K, Y] >= 0)                       # curtailable allocation above a
    @variable(m, q[K, Y, T] >= 0)                    # curtailable demand served in hour t
    @variable(m, z[K, Y], Bin)                       # admitted

    m[:lds] = lds
    loads_at = [Int[] for _ in 1:net.nb]; for k in K; push!(loads_at[lds[k].bus], k); end
    wd(n, yi, t) = isempty(loads_at[n]) ? AffExpr(0.0) : sum(a[k, yi] + q[k, yi, t] for k in loads_at[n])
    add_network!(m, net, sc, Y, wd; shed=shed)
    x, p, shed = m[:x], m[:p], m[:shed]

    @constraint(m, cap[k=K, yi=Y], a[k, yi] + f[k, yi] <= lds[k].L[yi] * z[k, yi])     # (6), named for the duals
    for k in K, yi in Y
        ld = lds[k]
        @constraint(m, f[k, yi] <= ld.phi * ld.L[yi] * z[k, yi])                         # (7)
        ld.pclr || @constraint(m, a[k, yi] >= ld.theta * z[k, yi])                       # (8)
        for t in T
            @constraint(m, q[k, yi, t] <= f[k, yi])                                      # (9)
            sc.stress[t] && fix(q[k, yi, t], 0.0; force=true)                            # (9) off in stress hours
        end
        @constraint(m, sum(sc.Nt[t] * (f[k, yi] - q[k, yi, t]) for t in T if !sc.stress[t])
                       <= ld.eps * N_YR * (a[k, yi] + f[k, yi]))                         # (10)
    end
    for k in K, yi in 1:length(Y)-1                                                       # (11), consecutive years suffice
        @constraint(m, z[k, yi] <= z[k, yi+1])
        @constraint(m, a[k, yi] <= a[k, yi+1])
        @constraint(m, a[k, yi] + f[k, yi] <= a[k, yi+1] + f[k, yi+1])
    end

    if variant in (:i, :ii)
        for l in 1:net.nl, yi in Y; fix(x[l, yi], 0.0; force=true); end
    end
    if variant == :i
        for k in K, yi in Y
            fix(f[k, yi], 0.0; force=true)
            for t in T; fix(q[k, yi, t], 0.0; force=true); end
        end
    end

    # (1) welfare; shed penalty added so a background-demand shortfall shows up instead of infeasibility
    yr = sc.years
    @expression(m, served[yi=Y], sum(lds[k].v * (a[k, yi] + BETA * sum(sc.Nt[t] / N_YR * q[k, yi, t] for t in T)) for k in K))
    @expression(m, gencost[yi=Y], sum(sc.Nt[t] * net.gamma[g] * p[g, yi, t] for g in 1:net.ng, t in T))
    @expression(m, upgcost[yi=Y], sum(net.C[l] * x[l, yi] for l in 1:net.nl if net.expand[l]))
    @expression(m, shedcost[yi=Y], VOLL * sum(sc.Nt[t] * shed[n, yi, t] for n in 1:net.nb, t in T))
    @objective(m, Max, sum(DELTA[yr[yi]] * (served[yi] - gencost[yi] - upgcost[yi] - shedcost[yi]) for yi in Y))
    return m
end

# Fix z (and optionally a, f) and drop integrality; used for SQ welfare and for the request-cap duals
function fix_allocation!(m; z=nothing, a=nothing, f=nothing)
    for k in axes(m[:a], 1), yi in axes(m[:a], 2)
        if z !== nothing
            is_binary(m[:z][k, yi]) && unset_binary(m[:z][k, yi])
            fix(m[:z][k, yi], z[k, yi]; force=true)
        end
        a === nothing || fix(m[:a][k, yi], a[k, yi]; force=true)
        f === nothing || fix(m[:f][k, yi], f[k, yi]; force=true)
    end
    return m
end
