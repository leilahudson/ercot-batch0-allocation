# bz.jl — ERCOT's Batch Zero rule (BZ) as we model it
# ERCOT allocates peak demand with a SCOPF whose objective is not published. We use OPT (i)'s objective
# with the threshold (8) removed and z = 1, one LP per year: firm MW only, f = 0, x = 0, N-1 secure in
# every modeled hour. With equal v_k the LP usually has many optimal allocations, so each draw values
# load k at v_k (1 + w_k), w_k uniform on [0, BZ_TIE], to pick one. Then allocations below the threshold
# go to zero with no reallocation, and a_{k,y} <- min_{y' >= y} a_{k,y'}. Zeroing below a threshold and
# taking the minimum commute, so the order of the two steps does not matter. The result is each load's
# LPC. A PCLR is studied at its requested peak, so it also gets the rest, f = L - a, as curtailable MW
# that the dispatch serves as far as the network allows outside stress hours.

function run_bz(net::Network, sc::Scenario, lds::Vector{LargeLoad}; draw::Int=1, tag="BZ")
    K, nY = length(lds), length(sc.years)
    w = BZ_TIE .* rand(MersenneTwister(SEED + 1000 + draw), K)
    a_raw = zeros(K, nY)              # after the threshold, before the min-over-later-years step
    amax  = zeros(K, nY)              # SCOPF allocation, before the threshold
    for yi in 1:nY
        m = new_model()
        @variable(m, 0 <= a[k=1:K] <= lds[k].L[yi])
        loads_at = [Int[] for _ in 1:net.nb]; for k in 1:K; push!(loads_at[lds[k].bus], k); end
        wd(n, y, t) = isempty(loads_at[n]) ? AffExpr(0.0) : sum(a[k] for k in loads_at[n])
        add_network!(m, net, sc, yi:yi, wd)
        for l in 1:net.nl; fix(m[:x][l, yi], 0.0; force=true); end
        p, shed = m[:p], m[:shed]
        @objective(m, Max, sum(lds[k].v * (1 + w[k]) * a[k] for k in 1:K)
                           - sum(sc.Nt[t] * net.gamma[g] * p[g, yi, t] for g in 1:net.ng, t in 1:sc.T)
                           - VOLL * sum(sc.Nt[t] * shed[n, yi, t] for n in 1:net.nb, t in 1:sc.T))
        st = solve_n1!(m, net, sc, yi:yi; tag="$tag y$(sc.years[yi])")
        st == MOI.OPTIMAL || error("[$tag] SCOPF status $st")
        for k in 1:K
            amax[k, yi] = clamp(value(a[k]), 0.0, lds[k].L[yi])
            keep = lds[k].pclr ? amax[k, yi] > 1e-3 : amax[k, yi] >= lds[k].theta - 1e-6
            a_raw[k, yi] = keep ? amax[k, yi] : 0.0
        end
    end
    a_bz = min_later_years(a_raw)
    f_bz = [lds[k].pclr ? lds[k].L[yi] - a_bz[k, yi] : 0.0 for k in 1:K, yi in 1:nY]
    return a_bz, f_bz, amax
end

# a shortfall in a later year lowers the earlier ones too, as ERCOT does in Batch Zero (keeps (10))
function min_later_years(a_raw)
    a = copy(a_raw)
    for k in axes(a_raw, 1), yi in axes(a_raw, 2)
        a[k, yi] = minimum(a_raw[k, yi:end])
    end
    return a
end
