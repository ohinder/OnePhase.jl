using advanced_timer, JLD

MUMPS = false
include("parameters.jl")
include("utils/utils.jl")
include("linear_system_solvers/linear_system_solvers.jl")
include("kkt_system_solver/include.jl")
include("line_search.jl")
include("IPM_tools.jl")
include("toy_examples.jl")
include("one_phase.jl")
include("init.jl")
include("delta_strategy.jl")

my_par = Class_parameters()
