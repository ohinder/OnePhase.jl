include("../include.jl")

using JuMP, Ipopt

m = Model()
@variable(m, x >= 0)
@variable(m, y >= 0)
@objective(m, Min, x-y)
@NLconstraint(m, x^2 + y^2 == 1)

it = one_phase_solve(m);
@show it.point.x

IpoptSolve(nlp_raw)
