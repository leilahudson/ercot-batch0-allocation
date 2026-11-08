# ERCOT Batch Zero allocation

15.S09 project, MIT, Fall 2026.

Compares a queue-order study rule (SQ) with a coordinated model (CO) for allocating transmission
headroom to large loads in ERCOT's Batch Zero. CO picks firm MW, curtailable MW and line upgrades
for all loads at once, N-1 secure. Three versions: (i) firm only, (ii) + curtailable, (iii) + upgrades.

## Files

- `config.jl` – settings, each can be overridden with an env var of the same name
- `batch0.jl` – loads everything, included by the run scripts
- `network.jl` – reads the case, builds PTDF and LODF
- `scenario.jl` – hours, demand, wind/solar, large-load requests
- `model.jl` – CO model and N-1 constraint generation
- `sq.jl` – status-quo rule
- `metrics.jl` – welfare, stranded MW, output tables
- `run_toy.jl` – 3-bus example
- `run_texas.jl` – Texas study
- `run.sbatch` – SLURM script for ORCD
- `prep/` – builds `data/*.csv` from the raw data
- `results/` – toy and Texas results

## Setup

Julia 1.12. Gurobi for the Texas runs, HiGHS is enough for the toy.

```
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Toy example

```
julia --project=. run_toy.jl
```

Output goes to `results/toy/`.

## Texas runs

Stride 2 (every second hour plus all stress hours), seed 1. Sites are at random buses in their
zone until bus coordinates are added.

```
R="--partition=mit_normal --cpus-per-task=8 --mem=96G --time=08:00:00"
MODELS=sq  HOUR_STRIDE=2 RUN_TAG=full_s2_sq sbatch $R run.sbatch
MODELS=co1 HOUR_STRIDE=2 RUN_TAG=full_s2_co sbatch $R run.sbatch
MODELS=co2 HOUR_STRIDE=2 RUN_TAG=s2_co2 WARM_FROM=results/full_s2_co/alloc_CO_i.csv sbatch $R run.sbatch
MODELS=co3 HOUR_STRIDE=2 RUN_TAG=s2_co3 WARM_FROM=results/s2_co2/alloc_CO_ii.csv sbatch $R run.sbatch
```

Each run writes `summary.csv`, `alloc_<model>.csv` and `loads.csv` to `results/<RUN_TAG>/`, plus
`bounds.csv` for CO (LP bound, final value, gap).

## Data

- `case_ACTIVSg2000.m` – Texas A&M synthetic grid, https://electricgrids.engr.tamu.edu
- zone load, wind and solar (2025, EIA-930) and data center sites (Baxtel, Cleanview) – from
  Thomas Lee's data bundle

Raw data isn't included. `prep/prepare_ercot_inputs.py` rebuilds the processed files from it.
