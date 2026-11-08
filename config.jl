# config.jl — run settings for the Batch Zero allocation study

_env(k, d::String)  = get(ENV, k, d)
_env(k, d::Int)     = parse(Int, get(ENV, k, string(d)))
_env(k, d::Float64) = parse(Float64, get(ENV, k, string(d)))
_env(k, d::Bool)    = lowercase(get(ENV, k, string(d))) == "true"

# solver
const SOLVER      = _env("SOLVER", "highs")        # "gurobi" on ORCD, "highs" for local tests
const THREADS     = _env("THREADS", 8)
const MIP_GAP     = _env("MIP_GAP", 5e-3)          # relative to an objective dominated by background generation cost
const TIME_LIMIT  = _env("TIME_LIMIT", 14400.0)    # s per solve; CO (iii) on 50 loads x 72 h needed more than 1 h
const MILP_TIME   = _env("MILP_TIME", 10800.0)
const START_TIME  = _env("START_TIME", 1200.0)      # s per solve when testing a start candidate; an infeasible candidate
                                                    # (too many loads forced to their threshold) is then dropped quickly       # s for CO's branch and bound; it starts from a feasible solution and
                                                    # bounds.csv records the gap to the LP bound, so a short run is still certified
const VERBOSE     = _env("VERBOSE", false)          # solver log on/off
const GRB_NODEFILE = _env("GRB_NODEFILE", 8.0)      # GB of branch-and-bound nodes held in memory before Gurobi writes them to disk (0 = never)
const GRB_METHOD  = _env("GRB_METHOD", -1)          # Gurobi Method for CO: -1 auto, 1 dual simplex, 2 barrier
const SQ_METHOD   = _env("SQ_METHOD", 1)            # SQ: dual simplex, warm-starts each load and each cut round

# data
const DATA_DIR    = joinpath(@__DIR__, "data")
const CASE_FILE   = _env("CASE_FILE", joinpath(DATA_DIR, "case_ACTIVSg2000.m"))
const PROFILE_CSV = _env("PROFILE_CSV", joinpath(DATA_DIR, "profiles.csv"))              # 3 real 2025 ERCOT days, from prep/prepare_ercot_inputs.py
const ZONE_CSV    = _env("ZONE_CSV", joinpath(DATA_DIR, "zone_requests.csv"))            # zonal shares, only for LOADS_SOURCE=sample
const DC_CSV      = _env("DC_CSV", joinpath(DATA_DIR, "dc_sites.csv"))                   # ERCOT data center sites (Thomas Lee's data_bundle)
const ZPEAK_CSV   = _env("ZPEAK_CSV", joinpath(DATA_DIR, "zone_peaks_2025.csv"))         # 2025 weather-zone peaks
const COORD_CSV   = _env("COORD_CSV", joinpath(DATA_DIR, "activsg2000_bus_coords.csv"))  # bus_id,lat,lon; optional
const OUT_DIR     = _env("OUT_DIR", joinpath(@__DIR__, "results"))

# study years and hours
const YEARS       = [2028, 2030, 2032]              # \mathcal Y
const BASE_YEAR   = 2025                            # year of the load data; the case's Pd is taken as each zone's 2025 peak
const GROWTH      = _env("GROWTH", 0.01)            # background demand growth per yr
const HOUR_STRIDE = _env("HOUR_STRIDE", 1)          # keep every k-th hour of each day plus all stress hours; weights reassigned
const STRESS_FRAC = _env("STRESS_FRAC", 0.10)       # \mathcal S = top 10% of net load hours
const HYDRO_AF    = _env("HYDRO_AF", 0.5)           # hydro availability (no profile in the case)
const N_YR        = 8760.0                          # N^yr

# network
const ALPHA       = _env("ALPHA", 1.15)             # \alpha, emergency/normal rating (case has no rateB)
const MIN_KV_LOAD = _env("MIN_KV_LOAD", 115.0)      # large loads sit at buses at or above this kV
const UPG_COST    = _env("UPG_COST", 20_000.0)      # C_l for a median-length line [$/MW-yr], scaled by x_l
const UPG_MAX     = _env("UPG_MAX", 1.0)            # \bar x_l as a share of the normal rating

# large loads
const LOADS_SOURCE = _env("LOADS_SOURCE", "dc")     # "dc": real ERCOT data center sites; "sample": random requests
const N_LOADS     = _env("N_LOADS", 0)              # dc: 0 = every site, n = random subset of n sites (quick tests)
const MIN_SITE_MW = _env("MIN_SITE_MW", 75.0)       # dc: sites with at least this planned MW are studied loads
# dc: put existing data centers at their own buses (true) or leave them inside zone load (false).
# true needs real coordinates (COORD_CSV); at random buses, sites of up to 1 GW land on weak buses and
# the base case cannot serve background demand.
const EXISTING_AT_SITES = _env("EXISTING_AT_SITES", false)
const TOTAL_GW    = _env("TOTAL_GW", 40.0)          # sample only: total studied request in 2032
const P_PCLR      = _env("P_PCLR", 0.3)             # share of loads electing PCLR
const SEED        = _env("SEED", 1)
const THRESH_MW   = 100.0                           # \theta_k rule (ERCOT comments, May 2026)
const THRESH_SHR  = 0.9

# status-quo study rule (sq.jl)
# "prop": study practice. Secure least-cost base dispatch without the studied loads; each load is
#   balanced by raising all generators in proportion to their headroom; a load is limited only by
#   (post-)contingency constraints where its distribution factor is at least SQ_DFAX_MIN.
# "lp":   each load gets the most MW over every possible dispatch (the original rule). One load can
#   use up all redispatch flexibility and block every later load, so use only for comparison.
const SQ_RULE     = _env("SQ_RULE", "prop")
const SQ_DFAX_MIN = _env("SQ_DFAX_MIN", 0.03)       # materiality cutoff on MW flow per MW of load

# economics (objective)
const V_FIRM      = _env("V_FIRM", 1.0e6)           # v_k [$/MW-yr], equal for all loads
const BETA        = _env("BETA", 0.7)               # \beta
const EPS_CURT    = _env("EPS_CURT", 0.0025)        # \varepsilon_k, Norris et al. 0.25%
const DELTA       = Dict(2028 => 1.0, 2030 => 1.0, 2032 => 1.0)   # \delta_y
# Background shed is a slack, NOT in the proposal; it should be 0 in every run. The penalties must be
# large: one MW shed at the right bus can relieve a binding line enough to admit 100+ MW of load,
# so a weak penalty lets the model trade background demand for large load.
const VOLL        = _env("VOLL", 1.0e6)             # CO [$/MWh]
const CO_SHED     = _env("CO_SHED", false)          # shed slack in CO; off saves ~40% of CO's variables (shed was 0 in every CO run)
const SQ_SHED_PEN = _env("SQ_SHED_PEN", 1.0e6)      # SQ, per avg MW shed (objective is MW of the load being studied)

# N-1 constraint generation 
const N1_TOL      = _env("N1_TOL", 0.01)            # MW, violation tolerance (keep small, see sq.jl)
const N1_MAXCUTS  = _env("N1_MAXCUTS", 20000)       # cuts added per round (5000 was hit in the first rounds on the full case)
const N1_ROUNDS   = _env("N1_ROUNDS", 30)
const N1_ON       = _env("N1_ON", true)             # false = intact limits only (debugging)
const N1_LOG      = _env("N1_LOG", false)           # print every cut round
# MILP stage: N-1 cuts as lazy constraints inside one branch-and-bound (Gurobi, GLPK). Without it each
# MILP round is a full re-solve, and on the full case each one ran into TIME_LIMIT.
const N1_LAZY     = _env("N1_LAZY", true)
