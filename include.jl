using advanced_timer, JLD

MUMPS = false
include("IPM/parameters.jl")
include("utils/utils.jl")
include("linear_system_solvers/linear_system_solvers.jl")
include("kkt_system_solver/include.jl")
include("line_search/line_search.jl")
include("IPM/IPM_tools.jl")
include("IPM/display_progress.jl")
include("IPM/one_phase.jl")
include("IPM/init.jl")
include("IPM/delta_strategy.jl")

my_par = Class_parameters()
