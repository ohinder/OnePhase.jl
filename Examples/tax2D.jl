# var c{i in A, j in B} >= clb[i,j];  # consumption of tax payer (i,j)
# var y{i in A, j in B} >= ylb[i,j];  # income      of tax payer (i,j)

#=
maximize f:
   sum{i in A, j in B}
      (lambda[i,j] * (log(c[i,j]) - (y[i,j]/w[i])^mu1[j] / mu1[j])
       - 0.5 * gamma * (c[i,j]^2 + y[i,j]^2));
=#
#=
subject to
   Incentive{(i,j) in T, (p,q) in T: if i=p then j!=q}:
     (log(c[i,j]) - (y[i,j]/w[i])^mu1[j] / mu1[j])
   - (log(c[p,q]) - (y[p,q]/w[i])^mu1[j] / mu1[j]) >= 0;
=#

#=
Technology:
   sum{i in A, j in B} lambda[i,j]*(y[i,j] - c[i,j]) >= 0;
=#
include("include.jl")

using JuMP
using Ipopt


srand(1)
n = 20

A = 1:n;
B = 1:n;
T = [(i,j) for i in A for j in B];
clb_ij = 1.0;
ylb_ij = 1.0;

lambda = rand(n,n);
w = rand(n) + 1.0;
mu = rand(n) + 1.0;
gamma = 1.0

m = Model(solver = IpoptSolver())
@variable(m, c[i in A,j in B] >= clb_ij)# for (i,j) in T)
@variable(m, y[i in A,j in B] >= ylb_ij)# for (i,j) in T )

@NLobjective(m, Min, -sum( lambda[i,j] * (log(c[i,j]) - (y[i,j]/w[i])^mu[j] / mu[j])  - 0.5 * gamma * (c[i,j]^2 + y[i,j]^2) for (i,j) = T) )

for (i,j) in T
  for (p,q) in T
    if i != p || j != p
      @NLconstraint(m,(log(c[i,j]) - (y[i,j]/w[i])^mu[j] / mu[j]) - (log(c[p,q]) - (y[p,q]/w[i])^mu[j] / mu[j]) >= 0)
    end
  end
end

#status = solve(m)

#println("Objective value: ", getobjectivevalue(m))
#println("x = ", getvalue(x))
#println("y = ", getvalue(y))



nlp_raw = MathProgNLPModel(m);
#status = solve(m)

#getvalue(x), getvalue(y)
#AbstractNLPModel
nlp = Class_CUTEst(nlp_raw);

timer = class_advanced_timer()
start_advanced_timer(timer)
#include("include.jl")
#intial_it = initial_point_satisfy_bounds(nlp, my_par)
start_advanced_timer(timer, "INIT")
my_par = Class_parameters()
intial_it = init(nlp, my_par, timer);
pause_advanced_timer(timer, "INIT")

pause_advanced_timer(timer)
print_timer_stats(timer)

start_advanced_timer(timer)


#intial_it = initial_point_generic(nlp, my_par, nlp_raw.meta.x0)

@assert(is_feasible(intial_it, my_par.comp_feas))
iter, status, hist, t, err = one_phase_IPM(intial_it, my_par, timer);

pause_advanced_timer(timer)

print_timer_stats(timer)

#finalize(nlp_raw)
