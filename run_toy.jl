# run_toy.jl — 3-bus toy example for Sec. 5: queue order strands headroom that CO uses
#
#        G1 (cheap, bus 1)
#        /            \
#   line 1-2        line 1-3        all lines x = 0.1, rating 200 MW, emergency 230 MW
#      /                \
#   bus 2 ---line 2-3--- bus 3
#   loads D (PCLR), A     loads B, C, small peaker G2 (50 MW)
#
# With dispatch fixed before an outage, losing either line out of bus 1 puts all of bus 1's export
# on the other, so buses 2 and 3 get at most 230 + 50 = 280 MW. Background uses 80 at peak, 40 off-peak.
# SQ, order D, A, B, C: D takes 120 MW firm. Each MW is met by G1 and G2 in proportion to headroom,
#   so 920/970 of it crosses a line out of bus 1, and A, B, C are each offered 38.2 MW, all below the
#   100 MW threshold, so all three get 0 (38.2 MW stranded, 41.8 MW lost to the fixed response).
# CO (i) fills the 200 MW of peak headroom (e.g. D 80 + B 120; which loads is not unique).
# CO (ii) adds 40 MW of curtailable load off-peak (D's PCLR share).
# CO (iii) buys rating on lines 1-2 and 1-3; how much depends on C_l / v_k, see the sweep at the end.

include("batch0.jl")

net = make_network(; kv=fill(230.0, 3), area=[1, 1, 1], Pd=[0.0, 40.0, 40.0],
                     fr=[1, 1, 2], to=[2, 3, 3], x=[0.1, 0.1, 0.1], rate=[200.0, 200.0, 200.0],
                     gbus=[1, 3], Pmax=[1000.0, 50.0], gamma=[20.0, 40.0], fuel=["ng", "ng"])

# two hours: a peak hour (stress, 10% of the year) and an off-peak hour
lf = [1.0, 0.5]; Nt = [876.0, 7884.0]
D  = zeros(net.nb, length(YEARS), 2)
for yi in 1:length(YEARS), t in 1:2; D[:, yi, t] .= net.Pd .* lf[t]; end
sc = Scenario(YEARS, 2, Nt, [true, false], D, repeat(net.Pmax, 1, 2), ["peak", "offpeak"], [1, 2])

#            k  bus  area  L by year             Lmax   pclr   phi  theta  eps       v
lds = [LargeLoad(1, 2, 1, [120.0, 120.0, 120.0], 120.0, true,  0.5, 100.0, EPS_CURT, V_FIRM),   # D, PCLR
       LargeLoad(2, 2, 1, [150.0, 150.0, 150.0], 150.0, false, 0.0, 100.0, EPS_CURT, V_FIRM),   # A
       LargeLoad(3, 3, 1, [120.0, 120.0, 120.0], 120.0, false, 0.0, 100.0, EPS_CURT, V_FIRM),   # B
       LargeLoad(4, 3, 1, [120.0, 120.0, 120.0], 120.0, false, 0.0, 100.0, EPS_CURT, V_FIRM)]   # C

rows = NamedTuple[]
a_sq, a_raw, amax = run_sq(net, sc, lds, [1, 2, 3, 4]; tag="SQ")
sol = evaluate_fixed(net, sc, lds, a_sq; tag="SQ eval")
strand = stranded(net, sc, lds, a_sq, a_raw, amax, [1, 2, 3, 4])
append!(rows, summary_rows("SQ", sol, sc, lds; amax=amax, stranded=strand))
println("SQ offered (before threshold), 2028: ", round.(amax[:, 1], digits=1), "  kept: ", round.(a_sq[:, 1], digits=1),
        "  stranded: ", round(strand[1], digits=1), " MW")

for var in (:i, :ii, :iii)
    m = build_co(net, sc, lds; variant=var)
    solve_n1!(m, net, sc, 1:length(YEARS); tag="CO-$var", relax_first=true)
    s = read_solution(m, net, sc, lds)
    append!(rows, summary_rows("CO_$var", s, sc, lds))
    println("CO-$var 2028  a = ", round.(s.a[:, 1], digits=1), "  f = ", round.(s.f[:, 1], digits=1),
            "  upgrades = ", round.(s.x[:, 1], digits=1))
    if var == :iii
        du = cap_duals(m, net, sc, lds, s.z)              # $/MW-yr of extra request, z fixed at the optimum
        # requests are flat across years and (11) ties the years together, so the LP dual can sit in any
        # year (degenerate); only the sum over years is unique: one more MW of request in every year
        println("CO-iii request-cap duals summed over years [\$/MW]: ", round.(vec(sum(du, dims=2)), digits=0))
    end
end

df = DataFrame(rows)
show(df[df.year .== YEARS[1], [:model, :firm_GW, :total_GW, :MW_equiv_GW, :n_admitted, :n_zeroed_by_threshold,
                                :stranded_GW, :cost_M, :welfare_M, :upg_GW, :shed_MWh]], allrows=true, allcols=true)
println()
mkpath(joinpath(OUT_DIR, "toy")); CSV.write(joinpath(OUT_DIR, "toy", "summary.csv"), df)

# ---- CO (iii) against the upgrade cost C_l (same for all three lines here) ----
# TOY_C_SWEEP="2e4,2e5,5e5,1e6,2e6" in $/MW-yr; v_k = V_FIRM
sweep = parse.(Float64, split(get(ENV, "TOY_C_SWEEP", "2e4,2e5,3e5,4e5,5e5,1e6"), ','))
srows = NamedTuple[]
C0 = copy(net.C)
for Cv in sweep
    net.C .= Cv
    m = build_co(net, sc, lds; variant=:iii)
    solve_n1!(m, net, sc, 1:length(YEARS); tag="sweep", relax_first=true)
    s = read_solution(m, net, sc, lds)
    r = summary_rows("CO_iii", s, sc, lds)[1]
    push!(srows, (C_per_MWyr=Cv, C_over_v=Cv / V_FIRM, upg_MW=r.upg_GW * 1e3, firm_MW=r.firm_GW * 1e3,
                  total_MW=r.total_GW * 1e3, n_admitted=r.n_admitted, welfare_M=r.welfare_M))
end
net.C .= C0
println("\nCO (iii) vs upgrade cost, 2028:")
show(DataFrame(srows), allrows=true, allcols=true); println()
CSV.write(joinpath(OUT_DIR, "toy", "upgrade_cost_sweep.csv"), DataFrame(srows))
