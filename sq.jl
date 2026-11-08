# sq.jl — status-quo queue-order rule (Appendix A)
# Each year on its own: loads in queue order, each gets the most firm MW (up to its request) that
# keeps the study's limits satisfied given the loads ahead of it, f = 0, x = 0. Below threshold -> 0
# (PCLRs are exempt). Afterwards a_{k,y} <- min_{y' >= y} a_{k,y'} (Planning Guide 9.3.2(6)).
# Two versions of "keeps the limits satisfied", chosen by SQ_RULE (config.jl):
#   run_sq_prop  fixed base dispatch + proportional generation response (default, study practice)
#   run_sq_lp    best possible dispatch for each load (original rule)

function queue_order(lds::Vector{LargeLoad}, rule::Symbol; seed=SEED)
    rule == :random   && return shuffle(MersenneTwister(seed + 100), collect(1:length(lds)))
    rule == :smallest && return sortperm([l.Lmax for l in lds])
    rule == :largest  && return sortperm([l.Lmax for l in lds], rev=true)
    error("unknown queue order $rule")
end

run_sq(net, sc, lds, order; tag="SQ") = SQ_RULE == "lp" ? run_sq_lp(net, sc, lds, order; tag=tag) :
                                                          run_sq_prop(net, sc, lds, order; tag=tag)

function run_sq_lp(net::Network, sc::Scenario, lds::Vector{LargeLoad}, order::Vector{Int}; tag="SQ")
    K, nY = length(lds), length(sc.years)
    a_raw = zeros(K, nY)              # allocation before the min-over-later-years step
    amax  = zeros(K, nY)              # most MW the LP offered, before the threshold
    for yi in 1:nY
        t0 = time()
        m = new_model(; method=SQ_METHOD)   # only a[k]'s bounds and the cuts change between solves
        @variable(m, a[k=1:K] >= 0)
        for k in 1:K; fix(a[k], 0.0; force=true); end
        loads_at = [Int[] for _ in 1:net.nb]; for k in 1:K; push!(loads_at[lds[k].bus], k); end
        wd(n, y, t) = isempty(loads_at[n]) ? AffExpr(0.0) : sum(a[k] for k in loads_at[n])
        add_network!(m, net, sc, yi:yi, wd)
        for l in 1:net.nl; fix(m[:x][l, yi], 0.0; force=true); end
        shed_avg = @expression(m, sum(sc.Nt[t] / N_YR * m[:shed][n, yi, t] for n in 1:net.nb, t in 1:sc.T))

        for k in order
            ld = lds[k]
            unfix(a[k]); set_lower_bound(a[k], 0.0); set_upper_bound(a[k], ld.L[yi])
            @objective(m, Max, a[k] - SQ_SHED_PEN * shed_avg)
            st = solve_n1!(m, net, sc, yi:yi; tag="$tag y$(sc.years[yi]) k$k")
            amax[k, yi] = st == MOI.OPTIMAL ? clamp(value(a[k]), 0.0, ld.L[yi]) : 0.0
            sh = st == MOI.OPTIMAL ? value(shed_avg) : NaN
            sh > 1e-3 && @printf("[%s] WARNING y%d k%d: %.3f avg MW of background shed\n", tag, sc.years[yi], k, sh)
            keep = ld.pclr ? amax[k, yi] > 1e-3 : amax[k, yi] >= ld.theta - 1e-6
            a_raw[k, yi] = keep ? amax[k, yi] : 0.0
            fix(a[k], a_raw[k, yi]; force=true)
        end
        @printf("[%s] %d: %.1f GW allocated to %d loads (%d zeroed by threshold), %.0fs\n", tag, sc.years[yi],
                sum(a_raw[:, yi]) / 1e3, sum(a_raw[:, yi] .> 0),
                sum(!lds[k].pclr && 1e-3 < amax[k, yi] < lds[k].theta for k in 1:K), time() - t0)
    end
    # a shortfall in a later year lowers the earlier ones too (keeps (11))
    return min_later_years(a_raw), a_raw, amax
end

# a shortfall in a later year lowers the earlier ones too (keeps (11), Planning Guide 9.3.2(6))
function min_later_years(a_raw)
    a_sq = copy(a_raw)
    for k in axes(a_raw, 1), yi in axes(a_raw, 2)
        a_sq[k, yi] = minimum(a_raw[k, yi:end])
    end
    return a_sq
end

# Study-practice SQ; for each year:
#  1. base case: least-cost dispatch p0 of background demand only, intact and N-1 limits (the
#     secure starting case a planner would study from), flows F0
#  2. each MW of a studied load at bus b is balanced by raising generator g by pi_g,t, its share of
#     total headroom Pbar - p0 in hour t, so its flow factor is d[l,t] = sum_g pi_g,t PTDF[l,b_g] - PTDF[l,b]
#  3. loads in queue order: the most MW such that every intact and post-contingency flow with a
#     factor of at least SQ_DFAX_MIN stays within its limit, in every hour, and total generation
#     headroom is not exceeded. The flows then include that load before the next one is studied.
# Closed form per load, no LP. Ignoring factors below SQ_DFAX_MIN is the usual materiality rule; it
# means many small impacts can add up, which shows up as shed in the SQ welfare evaluation.

# base case for year yi: secure least-cost dispatch without studied loads
function sq_base_case(net::Network, sc::Scenario, yi; tag="SQ base")
    m = new_model()
    add_network!(m, net, sc, yi:yi, (n, y, t) -> AffExpr(0.0))
    for l in 1:net.nl; fix(m[:x][l, yi], 0.0; force=true); end
    @objective(m, Min, sum(sc.Nt[t] * (net.gamma[g] * m[:p][g, yi, t]) for g in 1:net.ng, t in 1:sc.T) +
                       VOLL * sum(sc.Nt[t] * m[:shed][n, yi, t] for n in 1:net.nb, t in 1:sc.T))
    st = solve_n1!(m, net, sc, yi:yi; tag=tag)
    st == MOI.OPTIMAL || error("[$tag] base case status $st")
    p0 = [value(m[:p][g, yi, t]) for g in 1:net.ng, t in 1:sc.T]
    F0 = [value(m[:F][l, yi, t]) for l in 1:net.nl, t in 1:sc.T]
    shed = maximum(sum(value(m[:shed][n, yi, t]) for n in 1:net.nb) for t in 1:sc.T)
    shed > 1e-3 && @printf("[%s] WARNING: base case sheds up to %.1f MW of background demand in an hour\n", tag, shed)
    return p0, F0
end

# most MW one load can add given current flows F (nl x T), its flow factors d (nl x T) and the
# remaining generation headroom R (T); the limits are N1_TOL looser, matching the cut tolerance
function sq_max_mw(net::Network, F, d, R, cap)
    best = cap
    T = size(F, 2)
    for t in 1:T
        best = min(best, R[t])
        @inbounds for l in 1:net.nl                      # intact
            dl = d[l, t]
            abs(dl) < SQ_DFAX_MIN && continue
            room = dl > 0 ? net.rate[l] + N1_TOL - F[l, t] : net.rate[l] + N1_TOL + F[l, t]
            best = min(best, max(room, 0.0) / abs(dl))
        end
        @inbounds for c in 1:net.nl                      # post-contingency, outage of c
            net.cont_ok[c] || continue
            Fc, dc = F[c, t], d[c, t]
            for l in 1:net.nl
                l == c && continue
                L = net.LODF[l, c]
                dp = d[l, t] + L * dc
                abs(dp) < SQ_DFAX_MIN && continue
                Fp = F[l, t] + L * Fc
                lim = ALPHA * net.rate[l] + N1_TOL
                room = dp > 0 ? lim - Fp : lim + Fp
                best = min(best, max(room, 0.0) / abs(dp))
            end
        end
    end
    return max(best, 0.0)
end

function run_sq_prop(net::Network, sc::Scenario, lds::Vector{LargeLoad}, order::Vector{Int}; tag="SQ")
    K, nY, T = length(lds), length(sc.years), sc.T
    a_raw = zeros(K, nY); amax = zeros(K, nY)
    for yi in 1:nY
        t0 = time()
        p0, F = sq_base_case(net, sc, yi; tag="$tag base y$(sc.years[yi])")
        head = max.(sc.Pbar .- p0, 0.0)                         # ng x T
        R = vec(sum(head, dims=1))                              # remaining headroom per hour
        G = net.PTDF[:, net.gbus] * (head ./ max.(R', 1e-9))    # nl x T, flow per MW of proportional response
        for k in order
            ld = lds[k]
            d = G .- net.PTDF[:, ld.bus]                        # nl x T, flow per MW of load k
            amax[k, yi] = sq_max_mw(net, F, d, R, ld.L[yi])
            keep = ld.pclr ? amax[k, yi] > 1e-3 : amax[k, yi] >= ld.theta - 1e-6
            a_raw[k, yi] = keep ? amax[k, yi] : 0.0
            F .+= a_raw[k, yi] .* d; R .-= a_raw[k, yi]           # the next load is studied with this one in
        end
        @printf("[%s] %d: %.1f GW allocated to %d loads (%d zeroed by threshold), %.0fs\n", tag, sc.years[yi],
                sum(a_raw[:, yi]) / 1e3, sum(a_raw[:, yi] .> 0),
                sum(!lds[k].pclr && 1e-3 < amax[k, yi] < lds[k].theta for k in 1:K), time() - t0)
        flush(stdout)
    end
    return min_later_years(a_raw), a_raw, amax
end

# Stranded MW under the study rule: after SQ, the threshold-zeroed loads are studied again in queue
# order with no threshold, each seeing the ones before it, and their MW are summed. This counts the
# headroom they could have used without double counting what they share.
function stranded_prop(net::Network, sc::Scenario, lds, a_raw, amax, order; tag="SQ stranded")
    out = zeros(length(sc.years))
    for yi in eachindex(sc.years)
        Z = [k for k in order if !lds[k].pclr && 1e-3 < amax[k, yi] < lds[k].theta]
        isempty(Z) && continue
        p0, F = sq_base_case(net, sc, yi; tag="$tag base y$(sc.years[yi])")
        head = max.(sc.Pbar .- p0, 0.0); R = vec(sum(head, dims=1))
        G = net.PTDF[:, net.gbus] * (head ./ max.(R', 1e-9))
        for k in order                                          # rebuild flows with the SQ allocation
            a_raw[k, yi] > 0 || continue
            F .+= a_raw[k, yi] .* (G .- net.PTDF[:, lds[k].bus]); R .-= a_raw[k, yi]
        end
        for k in Z
            d = G .- net.PTDF[:, lds[k].bus]
            v = sq_max_mw(net, F, d, R, lds[k].L[yi])
            out[yi] += v; F .+= v .* d; R .-= v
        end
    end
    return out
end
