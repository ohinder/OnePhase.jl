
include("../benchmark.jl")
using CUTEst, OnePhase
#nlp_raw2 = CUTEstModel("RK23")
#nlp_raw2 = CUTEstModel("FLOSP2TH")

my_pars = OnePhase.Class_parameters()
nlp_raw2 = CUTEstModel("CAMSHAPE")

#nlp_raw2 = CUTEstModel("NET4")

if false
  my_pars.term.max_it = 100
  new_pars = OnePhase.autotune(nlp_raw2, my_pars)
end
#my_pars.term.tol_opt = 1e-6

iter, status, hist, t, err = OnePhase.one_phase_solve(nlp_raw2,my_pars);
#get_fval(iter)
#eval_f(nlp,iter.point.x)
#obj(nlp_raw2,iter.point.x) - get_fval(iter)
if false
using Ipopt
sol = IpoptSolve(nlp_raw2);
end

finalize(nlp_raw2)
