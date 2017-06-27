include("include.jl")

# INFEASIBLE PROBLEMS
#nlp_raw = CUTEstModel("NCVXQP8")
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

#nlp_raw = CUTEstModel("DISC2")
#nlp_raw = CUTEstModel("TFI1")
#nlp_raw = CUTEstModel("GPP")
#nlp_raw = CUTEstModel("ANTWERP")
#nlp_raw = CUTEstModel("HYDCAR20")
#nlp_raw = CUTEstModel("STEENBRD")
#nlp_raw = CUTEstModel("STEENBRC")
#nlp_raw = CUTEstModel("EXPFITC")
#nlp_raw = CUTEstModel("LAUNCH")
#nlp_raw = CUTEstModel("TRIMLOSS")

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
#nlp_raw = CUTEstModel("HYDCAR6")

#HVYCRASH,
#DISCS, EQC, HIMMELBJ,  PFIT1, PFIT3, SSEBNLN
#nlp_raw = CUTEstModel("DISCS")
#nlp_raw = CUTEstModel("HYDCAR6")
#nlp_raw = CUTEstModel("DISC2")
#nlp_raw = CUTEstModel("SSEBNLN")
#nlp_raw = CUTEstModel("BATCH")
#nlp_raw = CUTEstModel("SAWPATH")

#nlp_raw = CUTEstModel("GPP")
#mean(abs(grad(nlp_raw, nlp_raw.meta.x0))) #, maximum(abs(jac(nlp_raw, nlp_raw.meta.x0)))
#maximum(abs(grad(nlp_raw, nlp_raw.meta.x0)))
#mean(abs(jac(nlp_raw, nlp_raw.meta.x0)))
#finalize(nlp_raw)

## HARD PROBLEMS
#nlp_raw = CUTEstModel("ACOPP57")
#nlp_raw = CUTEstModel("ACOPP300")
#nlp_raw = CUTEstModel("LEAKNET")
#nlp_raw = CUTEstModel("TFI1")
#nlp_raw = CUTEstModel("TRAINH") # >> 1000, STRUGGLING, LINEAR SOLVER IS NOT V. GOOD # IPOPT 58
nlp_raw = CUTEstModel("AVION2") # HARD and poorly conditioned
#nlp_raw = CUTEstModel("A4X12") # HARD and poorly conditioned
#nlp_raw = CUTEstModel("CRESC100") # >> 100. Infinities!
#nlp_raw = CUTEstModel("AVION2")
#nlp_raw = CUTEstModel("QPNSTAIR")
#nlp_raw = CUTEstModel("YORKNET")

if false
using Ipopt
mp = NLPModels.NLPtoMPB(nlp_raw, IpoptSolver(print_level=5, tol=1e-8))
MathProgBase.optimize!(mp)
x = MathProgBase.getsolution(mp)
solver = MathProgBase.getrawsolver(mp)
end

begin
nlp = Class_CUTEst(nlp_raw)
## FEASIBLE (probably)
timer = class_advanced_timer()
start_advanced_timer(timer)
#include("include.jl")
#intial_it = initial_point_satisfy_bounds(nlp, my_par)
start_advanced_timer(timer, "INIT")
intial_it = init(nlp, my_par, timer)
pause_advanced_timer(timer, "INIT")

#intial_it = initial_point_generic(nlp, my_par, nlp_raw.meta.x0)

@assert(is_feasible(intial_it, my_par.comp_feas))
iter, status, hist, t, err = one_phase_IPM(intial_it, my_par, timer);

pause_advanced_timer(timer)

print_timer_stats(timer)

finalize(nlp_raw)
end
#
# aggressive steps do max LP step

if false
include("include.jl")
x = nlp_raw.meta.x0;
m = nlp_raw.meta.ncon;
@time for i = 1:20 obj(nlp_raw, x) end;
@time for i = 1:20 grad(nlp_raw, x) end;
@time for i = 1:20 cons(nlp_raw, x) end;
@time for i = 1:20 J = jac(nlp_raw, x) end;
@time for i = 1:20 p = jtprod(nlp_raw, x, randn(m))  end;
@time for i = 1:20  jac_coord(nlp_raw, x) end;
nlp = Class_CUTEst(nlp_raw);
@time for i = 1:20 eval_jac(nlp, zeros(10000)) end;
end
