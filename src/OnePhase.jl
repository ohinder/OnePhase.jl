 __precompile__()

module OnePhase

export Class_parameters, write_pars, one_phase_solve, autotune
export pd, rd, Eval_NaN_error, Class_CUTEst, init # utils

using advanced_timer, NLPModels, Compat, JuMP, CUTEst

#BLAS.set_num_threads(1) # this is because the server i uses does not correctly measure threads
USE_MUMPS = false
include("parameters.jl")
include("utils/utils.jl")
include("linear_system_solvers/linear_system_solvers.jl")
include("kkt_system_solver/include.jl")
include("line_search/line_search.jl")
include("IPM/ipm.jl")
include("init/init.jl")
include("JuMPinterface.jl")

end
