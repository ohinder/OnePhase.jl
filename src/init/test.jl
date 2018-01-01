include("../benchmark.jl")
#=
nlp_raw2 = CUTEstModel("AGG")
nlp = OnePhase.Class_CUTEst(nlp_raw2)
using advanced_timer
timer = advanced_timer.class_advanced_timer()
advanced_timer.start_advanced_timer(timer)
par = OnePhase.Class_parameters()
mu_ratio = OnePhase.gertz_init(nlp, par, timer)
advanced_timer.pause_advanced_timer(timer)
=#
#iter.point

finalize(nlp_raw2)
