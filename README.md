# ERCOT Batch Zero allocation

15.S09 project, MIT, Fall 2026.

Compares a queue-order study rule (SQ) with a coordinated model (CO) for allocating transmission
headroom to large loads in ERCOT's Batch Zero. CO picks firm MW, curtailable MW and line upgrades
for all loads at once, N-1 secure. Three versions: (i) firm only, (ii) + curtailable, (iii) + upgrades.

## Files

- `config.jl` – settings, each can be overridden with an env var of the same name
- `batch0.jl` – loads everything, included by the run script
- `network.jl` – builds the network, PTDF and LODF
- `scenario.jl` – hours and large-load requests
- `model.jl` – CO model and N-1 constraint generation
- `sq.jl` – status-quo rule
- `metrics.jl` – welfare and output tables
- `run_toy.jl` – 3-bus example
- run.sbatch – SLURM script for the cluster

## Run

Julia 1.12. HiGHS is enough, Gurobi also works.

```
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. run_toy.jl
SOLVER=gurobi julia --project=. run_toy.jl    # with Gurobi
```

Output goes to `results/toy/summary.csv` (the proposal's results table).
