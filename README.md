# Allocating transmission headroom to large loads in ERCOT's Batch Zero

Code and data for my 15.S09 project (MIT, Fall 2026). It compares a status-quo queue-order study
rule (SQ) with a coordinated allocation model (CO) that chooses every load's firm MW, curtailable MW
and line upgrades at once, with every line secure after any single outage. CO is solved in three
nested variants: (i) firm MW only, (ii) plus curtailable MW for PCLRs, (iii) plus line upgrades.

## Files

| File | What it does |
|---|---|
| `config.jl` | Every setting (solver, data paths, parameters). Each can be overridden with an env var of the same name. |
| `batch0.jl` | Loads config, the solver and the model files. Included by the run scripts. |
| `network.jl` | Reads the MATPOWER case, builds PTDF and LODF, makes the toy network. |
| `scenario.jl` | Hours, background demand, wind and solar availability, stress hours, large-load requests. |
| `model.jl` | Network constraints, N-1 constraint generation, CO (eq. 1–12), the MILP start heuristic and lazy N-1 cuts. |
| `sq.jl` | The status-quo rule (base case, proportional generation response, queue order, threshold). |
| `metrics.jl` | SQ welfare evaluation, stranded MW, request-cap duals, summary tables and CSVs. |
| `run_toy.jl` | The three-bus example (proposal Section 5, final report Appendix A). |
| `run_texas.jl` | The Texas study on ACTIVSg2000. |
| `run.sbatch` | SLURM script used on MIT ORCD. |
| `prep/` | Scripts that build `data/*.csv` from the raw data bundle. |
| `results/toy/` | Toy results (`summary.csv`, `upgrade_cost_sweep.csv`). |
| `results/texas_s2/` | Texas results behind the report's main table and figure, and the figure script. |

## Setup

Julia 1.12, JuMP, Gurobi (license needed for the Texas runs) and HiGHS (free, enough for the toy).

```
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Reproduce the toy example

```
julia --project=. run_toy.jl                    # HiGHS
SOLVER=gurobi julia --project=. run_toy.jl      # Gurobi
```

Writes `results/toy/summary.csv` (proposal Table 1) and `results/toy/upgrade_cost_sweep.csv`.

## Reproduce the Texas results

All runs use every second hour plus all stress hours (`HOUR_STRIDE=2`), seed 1, and sites at random
buses of at least 115 kV in their zone (no bus coordinates yet). On ORCD each was one job:

```
R="--partition=mit_normal --cpus-per-task=8 --mem=96G --time=08:00:00"
MODELS=sq  HOUR_STRIDE=2 RUN_TAG=full_s2_sq sbatch $R run.sbatch
MODELS=co1 HOUR_STRIDE=2 RUN_TAG=full_s2_co sbatch $R run.sbatch
MODELS=co2 HOUR_STRIDE=2 RUN_TAG=s2_co2 WARM_FROM=results/full_s2_co/alloc_CO_i.csv sbatch $R run.sbatch
MODELS=co3 HOUR_STRIDE=2 RUN_TAG=s2_co3 WARM_FROM=results/s2_co2/alloc_CO_ii.csv  sbatch $R run.sbatch
```

Each run writes `results/<RUN_TAG>/summary.csv` (one row per model and year), `alloc_<model>.csv`
(per load and year), `loads.csv`, and for CO `bounds.csv` (LP relaxation bound, start, final value,
Gurobi bound and gaps). `results/texas_s2/results_texas_s2.csv` collects the rows used in the
report, and `make_fig.py` draws Figure 1 from it.

## Data

| File | Source |
|---|---|
| `data/case_ACTIVSg2000.m` | Texas A&M synthetic grid ACTIVSg2000, https://electricgrids.engr.tamu.edu |
| `data/ercot_hourly_2025.csv`, `profiles.csv`, `zone_peaks_2025.csv` | 2025 EIA-930 zone load and ERCOT wind and solar capacity factors, via Thomas Lee's `data_bundle` |
| `data/dc_sites.csv` | ERCOT data center sites (Baxtel and Cleanview estimates), via Thomas Lee's `data_bundle` |

`prep/prepare_ercot_inputs.py` rebuilds the processed files from the raw bundle in `data/raw/`
(not included).

## How CO is solved

The N-1 limits are generated, not written out. Rounds on the LP relaxation add the worst violated
limit for each line, year and hour until none is violated, which gives an upper bound. A feasible
solution is then built by fixing admissions (the previous variant's, or the relaxation rounded) and
solving that LP with its own N-1 rounds. The MILP starts from it, with every integer solution
checked against all outages through lazy constraints. If branch and bound does not improve on the
start, the start is reported. `bounds.csv` records the certified gap.
