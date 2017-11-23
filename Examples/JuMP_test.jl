include("../include.jl")

using JuMP, Ipopt

m = Model()
@variable(m, x >= 0)
@variable(m, y >= 0)
@variable(m, z >= 0)

@objective(m, Min, x-y)
@NLconstraint(m, y <= 1)

@NLconstraint(m, x^2 + y^2 <= 1)
@NLconstraint(m, x^2 + y^2 >= 1)


nlp_raw = MathProgNLPModel(m);
my_pars = Class_parameters();
my_pars.term.tol_opt = 1e-30
it = one_phase_solve(nlp_raw,my_pars);
it.point.x
