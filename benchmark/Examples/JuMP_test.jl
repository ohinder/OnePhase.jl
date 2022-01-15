include("../benchmark.jl")

using JuMP #, Ipopt

m = Model(solver=OnePhaseSolver(init!mu_scale=0.01))
@variable(m, x >= 0)
@variable(m, y >= 0)
@objective(m, Min, x-y)
@NLconstraint(m, x^2 + y^2 <= 1)

solve(m)

getobjectivevalue(m)
