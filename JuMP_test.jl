include("include.jl")

using JuMP, Ipopt

m = Model(solver = IpoptSolver())
@variable(m, x >= 0)
@variable(m, y >= 0)
@objective(m, Max, 5x + 22y)
@NLconstraint(m, x + y^3 <= 1)

nlp_raw = MathProgNLPModel(m)
#status = solve(m)

#getvalue(x), getvalue(y)
#AbstractNLPModel
nlp = Class_CUTEst(nlp_raw)

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

#finalize(nlp_raw)
