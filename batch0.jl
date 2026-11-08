# batch0.jl — loads config, solver and model files

using Printf, Random
include("config.jl")

# Gurobi is loaded only when SOLVER=gurobi. These stay separate top-level lines so the
# newly loaded package is callable on the next line (Julia world age).
SOLVER == "gurobi" && @eval using Gurobi
const GRB_ENV = Ref{Any}(nothing)
SOLVER == "gurobi" && (GRB_ENV[] = Gurobi.Env())     # one license checkout for the whole run
SOLVER == "glpk" && @eval using GLPK                 # local tests of the lazy N-1 callback (HiGHS has none)

include("network.jl")
include("scenario.jl")
include("model.jl")
include("sq.jl")
include("metrics.jl")
