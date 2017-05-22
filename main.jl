include("include.jl")

# INFEASIBLE PROBLEMS
#nlp_raw = CUTEstModel("JUNKTURN")
#nlp_raw = CUTEstModel("DRCAVTY3") # seems to be feasible, IPOPT struggles
#nlp_raw = CUTEstModel("MODEL")

#nlp_raw = CUTEstModel("PT") # 13 ITS
#nlp_raw = CUTEstModel("AGG") # 153 ITS
#nlp_raw = CUTEstModel("KISSING") # 180 ITS
#nlp_raw = CUTEstModel("KISSING2") # 151 ITS
#nlp_raw = CUTEstModel("FLETCHCR")
#nlp_raw = CUTEstModel("GENHUMPS")
#nlp_raw = CUTEstModel("ZIGZAG") # 22 ITS
#nlp_raw = CUTEstModel("HVYCRASH") # 39 ITS
#nlp_raw = CUTEstModel("TRAINF") # 140 ITS or 81 ITS if AGG starts
#nlp_raw = CUTEstModel("ARTIF") # 14 ITS
#nlp_raw = CUTEstModel("AVGASB") # 9 ITS
#nlp_raw = CUTEstModel("HYDROELM")
#nlp_raw = CUTEstModel("DALLASM")
#nlp_raw = CUTEstModel("GPP")

#nlp_raw = CUTEstModel("GPP")


## HARD PROBLEMS
#nlp_raw = CUTEstModel("ACOPP57")
#nlp_raw = CUTEstModel("ACOPP300")

#nlp_raw = CUTEstModel("TRAINH") # >> 1000, STRUGGLING, LINEAR SOLVER IS NOT V. GOOD # IPOPT 58
#nlp_raw = CUTEstModel("AVION2") # HARD and poorly conditioned
#nlp_raw = CUTEstModel("A4X12") # HARD and poorly conditioned
nlp_raw = CUTEstModel("CRESC100") # >> 100. Infinities!
#nlp_raw = CUTEstModel("CHAIN")


if false
using Ipopt
mp = NLPModels.NLPtoMPB(nlp_raw, IpoptSolver(print_level=5,nlp_scaling_method="none"))
MathProgBase.optimize!(mp)
x = MathProgBase.getsolution(mp)
#MathProgBase.getdual(mp)
finalize(nlp_raw)
end

#nlp_raw = CUTEstModel("GENHUMPS")
nlp = Class_CUTEst(nlp_raw)
#my_par.tol = 1e-8

reset_advanced_timer()
start_advanced_timer()
#include("include.jl")
intial_it = initial_point_satisfy_bounds(nlp, my_par)
#intial_it = initial_point_generic(nlp, my_par)
iter, status, hist = one_phase_IPM(intial_it, my_par);
pause_advanced_timer()
print_timer_stats()

finalize(nlp_raw)
