# ERCOT Batch Zero allocation

15.S09 project, MIT, Fall 2026.

How should constrained transmission headroom be divided among competing large loads? This repo
compares three rules on the same network.

- **SQ**, the queue-order study rule ERCOT used before Batch Zero. Loads are studied one at a time,
  each against a fixed base dispatch, and an allocation below the 100 MW threshold is set to zero.
- **BZ**, ERCOT's Batch Zero procedure as we model it. A security-constrained OPF allocates MW to all
  loads at once, then allocations below the threshold are set to zero with no reallocation. ERCOT does
  not publish the SCOPF objective, so ties are broken with small random weights over `BZ_DRAWS` draws.
- **CO**, a coordinated model that chooses firm MW, curtailable MW (PCLRs only) and line upgrades for
  all loads at once, N-1 secure, with the threshold inside the optimization. Three versions, (i) firm
  only, (ii) adds curtailable MW, (iii) adds upgrades.

## Files

- `config.jl` – settings, each can be overridden with an environment variable of the same name
- `batch0.jl` – loads everything, included by the run script
- `network.jl` – MATPOWER reader, PTDF and LODF
- `scenario.jl` – representative hours and large-load requests
- `model.jl` – CO model and N-1 constraint generation
- `sq.jl` – SQ rule
- `bz.jl` – BZ rule
- `metrics.jl` – welfare and output tables
- `run_toy.jl` – the 3-bus example in the proposal

## Run

Julia 1.12. HiGHS is enough for the toy.

```
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. run_toy.jl
```

With Gurobi, install Gurobi.jl in your environment (it needs a license) and run with `SOLVER=gurobi`.

The toy writes `results/toy/summary.csv` (Table 1 of the proposal), `results/toy/bz_draws.csv`
(one row per BZ tie-breaking draw) and `results/toy/upgrade_cost_sweep.csv`.

The Texas study uses the ACTIVSg2000 case and data that are not in this repo (see `config.jl` for the
expected files under `data/`).
