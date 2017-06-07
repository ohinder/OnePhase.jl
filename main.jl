include("include.jl")

# INFEASIBLE PROBLEMS
#nlp_raw = CUTEstModel("JUNKTURN")
#nlp_raw = CUTEstModel("DRCAVTY3") # seems to be feasible, IPOPT struggles
#nlp_raw = CUTEstModel("MODEL")
#nlp_raw = CUTEstModel("KTMODEL")

#nlp_raw = CUTEstModel("QPCBOEI1")
#nlp_raw = CUTEstModel("PT") # 13 ITS
#nlp_raw = CUTEstModel("AGG") # 153 ITS
#nlp_raw = CUTEstModel("KISSING") # 180 ITS
#nlp_raw = CUTEstModel("KISSING2") # 151 ITS
#nlp_raw = CUTEstModel("FLETCHCR")
#nlp_raw = CUTEstModel("GENHUMPS")
nlp_raw = CUTEstModel("ZIGZAG") # 22 ITS
#nlp_raw = CUTEstModel("HVYCRASH") # 39 ITS
#nlp_raw = CUTEstModel("TRAINF") # 140 ITS or 81 ITS if AGG starts
#nlp_raw = CUTEstModel("ARTIF") # 14 ITS
#nlp_raw = CUTEstModel("AVGASB") # 9 ITS
#nlp_raw = CUTEstModel("HYDROELM")
#nlp_raw = CUTEstModel("AIRPORT")
#nlp_raw = CUTEstModel("QPCSTAIR")

#nlp_raw = CUTEstModel("DALLASM")
#nlp_raw = CUTEstModel("GPP")
#nlp_raw = CUTEstModel("ANTWERP")

#nlp_raw = CUTEstModel("STEENBRD")
#nlp_raw = CUTEstModel("STEENBRC")
#nlp_raw = CUTEstModel("EXPFITC")
#nlp_raw = CUTEstModel("LAUNCH")
#nlp_raw = CUTEstModel("HAIFAM")
#nlp_raw = CUTEstModel("LAKES")
#nlp_raw = CUTEstModel("ELATTAR")

#nlp_raw = CUTEstModel("CRESC50")




#nlp_raw = CUTEstModel("GPP")
#mean(abs(grad(nlp_raw, nlp_raw.meta.x0))) #, maximum(abs(jac(nlp_raw, nlp_raw.meta.x0)))
#maximum(abs(grad(nlp_raw, nlp_raw.meta.x0)))
#mean(abs(jac(nlp_raw, nlp_raw.meta.x0)))
#finalize(nlp_raw)

## HARD PROBLEMS
#nlp_raw = CUTEstModel("ACOPP57")
#nlp_raw = CUTEstModel("ACOPP300")

#nlp_raw = CUTEstModel("TRAINH") # >> 1000, STRUGGLING, LINEAR SOLVER IS NOT V. GOOD # IPOPT 58
#nlp_raw = CUTEstModel("AVION2") # HARD and poorly conditioned
#nlp_raw = CUTEstModel("A4X12") # HARD and poorly conditioned
#nlp_raw = CUTEstModel("CRESC100") # >> 100. Infinities!
#nlp_raw = CUTEstModel("CHAIN")
#nlp_raw = CUTEstModel("QPNSTAIR")


if false
using Ipopt
mp = NLPModels.NLPtoMPB(nlp_raw, IpoptSolver(print_level=5))
MathProgBase.optimize!(mp)
x = MathProgBase.getsolution(mp)
solver = MathProgBase.getrawsolver(mp)

#this_info = problem_summary()
#set_cutest_info_ipopt!(this_info, mp.inner, nlp_raw, x)
#finalize(nlp_raw)

#=
@show max(norm(mp.inner.mult_x_L,Inf), norm(mp.inner.mult_x_U,Inf), norm(mp.inner.mult_g,Inf))
num_vars = length(nlp_raw.meta.lvar)
num_cons = length(nlp_raw.meta.lcon)
#x_true = x[(end - num_vars):end]
x_true = x[1:num_vars]
#x_true = deepcopy(x)
obj(nlp_raw, x_true)
a = cons(nlp_raw, x_true);
primal_vio = max(0.0, maximum(nlp_raw.meta.lvar - x_true), maximum(x_true - nlp_raw.meta.uvar), maximum(nlp_raw.meta.lcon - a), maximum(a - nlp_raw.meta.ucon))
@show primal_vio
#@show max(0.0,)
dual = grad(nlp_raw, x_true) + jac(nlp_raw, x_true)' * mp.inner.mult_g + mp.inner.mult_x_U - mp.inner.mult_x_L;
@show norm(dual,Inf)=#
end
#eval_h_orginal = deepcopy(mp.inner.eval_h)

#nlp_raw = CUTEstModel("GENHUMPS")
nlp = Class_CUTEst(nlp_raw)
#my_par.tol = 1e-8

reset_advanced_timer()
start_advanced_timer()
#include("include.jl")
#intial_it = initial_point_satisfy_bounds(nlp, my_par)
intial_it = init(nlp, my_par)

if false
    convert_to_homog!(intial_it, my_par);
    convert_to_prox!(intial_it, my_par);
    intial_it.nlp.centre_point[end] = 1.0
    intial_it.nlp.lambda_vec[1:(end-1)] = 0.0
    intial_it.nlp.lambda_vec[end] = 10.0
end

@assert(is_feasible(intial_it, my_par.comp_feas))
iter, status, hist, t, err = one_phase_IPM(intial_it, my_par);
pause_advanced_timer()
print_timer_stats()

finalize(nlp_raw)
