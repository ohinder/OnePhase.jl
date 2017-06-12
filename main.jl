include("include.jl")

# INFEASIBLE PROBLEMS
#nlp_raw = CUTEstModel("JUNKTURN")
#nlp_raw = CUTEstModel("DRCAVTY3") # seems to be feasible, IPOPT struggles
#nlp_raw = CUTEstModel("MODEL")
#nlp_raw = CUTEstModel("FLOSP2HL")
#nlp_raw = CUTEstModel("WOODSNE")
#nlp_raw = CUTEstModel("FLOSP2HM")
#nlp_raw = CUTEstModel("CHNRSBNE")
#nlp_raw = CUTEstModel("KTMODEL")
#nlp_raw = CUTEstModel("CONT6-QQ")

#nlp_raw = CUTEstModel("QPCBOEI1")
#nlp_raw = CUTEstModel("PT") # 13 ITS
#nlp_raw = CUTEstModel("AGG") # 153 ITS
#nlp_raw = CUTEstModel("KISSING") # 180 ITS
#nlp_raw = CUTEstModel("KISSING2") # 151 ITS
#nlp_raw = CUTEstModel("FLETCHCR")
#nlp_raw = CUTEstModel("GENHUMPS")
#nlp_raw = CUTEstModel("ZIGZAG") # 22 ITS
#nlp_raw = CUTEstModel("HVYCRASH") # 39 ITS
#nlp_raw = CUTEstModel("TRAINF") # 140 ITS or 81 ITS if AGG starts
#nlp_raw = CUTEstModel("ARTIF") # 14 ITS, IPOPT infeasible
#nlp_raw = CUTEstModel("AVGASB") # 9 ITS
#nlp_raw = CUTEstModel("HYDROELM")
#nlp_raw = CUTEstModel("AIRPORT")
#nlp_raw = CUTEstModel("QPCSTAIR")

#nlp_raw = CUTEstModel("TFI1")
#nlp_raw = CUTEstModel("GPP")
#nlp_raw = CUTEstModel("ANTWERP")

#nlp_raw = CUTEstModel("STEENBRD")
#nlp_raw = CUTEstModel("STEENBRC")
#nlp_raw = CUTEstModel("EXPFITC")
#nlp_raw = CUTEstModel("LAUNCH")
#nlp_raw = CUTEstModel("HAIFAM")
#nlp_raw = CUTEstModel("ACOPR118")
#nlp_raw = CUTEstModel("LAKES")
#nlp_raw = CUTEstModel("ELATTAR")
#nlp_raw = CUTEstModel("ELEC")
#nlp_raw = CUTEstModel("CRESC50")
#nlp_raw = CUTEstModel("TRO4X4")
#nlp_raw = CUTEstModel("METHANL8")
#nlp_raw = CUTEstModel("GROUPING")
#nlp_raw = CUTEstModel("ARWHDNE")

#HVYCRASH,
#DISCS, EQC, HIMMELBJ,  PFIT1, PFIT3, SSEBNLN
#nlp_raw = CUTEstModel("SSEBNLN")


#nlp_raw = CUTEstModel("GPP")
#mean(abs(grad(nlp_raw, nlp_raw.meta.x0))) #, maximum(abs(jac(nlp_raw, nlp_raw.meta.x0)))
#maximum(abs(grad(nlp_raw, nlp_raw.meta.x0)))
#mean(abs(jac(nlp_raw, nlp_raw.meta.x0)))
#finalize(nlp_raw)

## HARD PROBLEMS
#nlp_raw = CUTEstModel("ACOPP57")
#nlp_raw = CUTEstModel("ACOPP300")
nlp_raw = CUTEstModel("LEAKNET")

#nlp_raw = CUTEstModel("TRAINH") # >> 1000, STRUGGLING, LINEAR SOLVER IS NOT V. GOOD # IPOPT 58
#nlp_raw = CUTEstModel("AVION2") # HARD and poorly conditioned
#nlp_raw = CUTEstModel("A4X12") # HARD and poorly conditioned
#nlp_raw = CUTEstModel("CRESC100") # >> 100. Infinities!
#nlp_raw = CUTEstModel("CHAIN")
#nlp_raw = CUTEstModel("QPNSTAIR")
#nlp_raw = CUTEstModel("YORKNET")

if true
using Ipopt
mp = NLPModels.NLPtoMPB(nlp_raw, IpoptSolver(print_level=5, tol=1e-8))
MathProgBase.optimize!(mp)
x = MathProgBase.getsolution(mp)
solver = MathProgBase.getrawsolver(mp)
end

nlp = Class_CUTEst(nlp_raw)

## FEASIBLE (probably)
if true
reset_advanced_timer()
start_advanced_timer()
#include("include.jl")
#intial_it = initial_point_satisfy_bounds(nlp, my_par)
intial_it = init(nlp, my_par)
@show norm(intial_it.point.x,2)
#intial_it = initial_point_generic(nlp, my_par, nlp_raw.meta.x0)

@assert(is_feasible(intial_it, my_par.comp_feas))
iter, status, hist, t, err = one_phase_IPM(intial_it, my_par);
pause_advanced_timer()
print_timer_stats()
end

## INFEASIBLE
if false
reset_advanced_timer()
start_advanced_timer()

c = ones(length(nlp_raw.meta.x0))
nlp_raw_infeas = Class_infeas_NLP(nlp_raw, -2.0, c)
nlp_infeas = Class_CUTEst(nlp_raw_infeas)

intial_it = init(nlp_infeas, my_par)

@assert(is_feasible(intial_it, my_par.comp_feas))


iter, status, hist, t, err = one_phase_IPM(intial_it, my_par);
pause_advanced_timer()
print_timer_stats()

if false
using Ipopt
mp = NLPModels.NLPtoMPB(nlp_raw_infeas, IpoptSolver(print_level=8))
MathProgBase.optimize!(mp)
x = MathProgBase.getsolution(mp)
solver = MathProgBase.getrawsolver(mp)
end
end




finalize(nlp_raw)
