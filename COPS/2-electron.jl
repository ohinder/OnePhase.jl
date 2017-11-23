include("../include.jl")

using JuMP #,

## OBJ DOES NOT MATCH



function create_2(n_p::Int64)
  m = Model()
  #@variable(m, x[i=1:n_p], start=(i - 0.5 * n_p)/n_p)
  #@variable(m, y[i=1:n_p], start=(i - 0.5 * n_p)/n_p)
  #@variable(m, z[i=1:n_p], start=(i - 0.5 * n_p)/n_p)
  @variable(m, x[i=1:n_p], start=rand())
  @variable(m, y[i=1:n_p], start=rand())
  @variable(m, z[i=1:n_p], start=rand())

  @NLobjective(m, Min, sum(sum(((x[i] - x[j])^2 + (y[i] - y[j])^2 + (z[i] - z[j])^2)^(-0.5) for j = (i+1):n_p) for i = 1:(n_p - 1)))

  for i = 1:n_p
    @NLconstraint(m, x[i]^2 + y[i]^2 + z[i]^2 == 1.0)
  end

  return m;
end

srand(1)
m = create_2(100)

nlp_raw = MathProgNLPModel(m);
my_pars = Class_parameters();
it = one_phase_solve(nlp_raw,my_pars);
it.point.x

#=
using Ipopt
setsolver(m,IpoptSolver())
solve(m)
=#
