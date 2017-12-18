include("../include.jl")

using JuMP, Ipopt

m = Model()
@variable(m, x[1:3])

eps = 1e-4
scaler = 1.1
@objective(m, Min, sum(x[i] for i = 1:3) )
@NLconstraint(m, scaler * x[1] <= 0)
@NLconstraint(m, x[1] >= eps)
@NLconstraint(m, scaler^2 * x[2] >= 0)
@NLconstraint(m, x[2] <= eps)
@NLconstraint(m, scaler^3 * x[3] >= 0)
@NLconstraint(m, x[3] <= eps)

nlp_raw = MathProgNLPModel(m);
my_pars = Class_parameters();
it = one_phase_solve(nlp_raw,my_pars);


IpoptSolve(nlp_raw)

#it.point.x
