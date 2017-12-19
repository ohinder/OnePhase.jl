include("../include.jl")

#nlp_raw = CUTEstModel("KISSING") # 180 ITS
#nlp_raw = CUTEstModel("KISSING2") # 151 ITS
nlp_raw = CUTEstModel("KISSING2","-param","m=100", "n=8")
#nlp_raw = CUTEstModel("CHANDHEU")
if true
using Ipopt
solver = IpoptSolver(print_level=5, max_iter=3000, bound_relax_factor=0.0, nlp_scaling_method="none")
#,mehrotra_algorithm="yes") #, tol_dual_abs=1e-6)
#solver = IpoptSolver(print_level=5, tol=1e-8)
mp = NLPModels.NLPtoMPB(nlp_raw, solver)
MathProgBase.optimize!(mp)
@show norm(mp.inner.mult_g, Inf)
#y = MathProgBase.getdual(mp)
solver = MathProgBase.getrawsolver(mp)
#finalize(nlp_raw)
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

@assert(is_feasible(intial_it, my_par.ls.comp_feas))
iter, status, hist, t, err = one_phase_IPM(intial_it, my_par, timer);

pause_advanced_timer(timer)

print_timer_stats(timer)

end
finalize(nlp_raw)

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
