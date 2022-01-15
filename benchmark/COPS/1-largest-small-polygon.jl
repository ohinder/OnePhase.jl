include("../include.jl")

using JuMP, Ipopt



function create_1(n_v::Int64)
  m = Model()
  @variable(m, 0.0 <= theta[i=1:n_v] <= pi)
  @variable(m, 0.0 <= r[i=1:n_v] <= 1.0)

  @NLobjective(m, Min, -0.5 * sum(r[i+1] * r[i] * sin(theta[i+1] - theta[i]) for i = 1:(n_v - 1)))
  for i = 1:(n_v-1)
    for j = (i+1):n_v
      @NLconstraint(m, r[i]^2 + r[j]^2 - 2 * r[i] * r[j] * cos(theta[i] - theta[j]) <= 1)
    end
    @NLconstraint(m, theta[i] <= theta[i+1])
  end
  @constraint(m, r[n_v] == 0.0)
  @constraint(m, theta[n_v] == pi)

  return m;
end

m = create_1(20)
nlp_raw = MathProgNLPModel(m);
my_pars = Class_parameters();
it = one_phase_solve(nlp_raw,my_pars);
it.point.x
#=
setsolver(m,IpoptSolver())
solve(m)
=#
