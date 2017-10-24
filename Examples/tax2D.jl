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
include("../include.jl")

using JuMP
using Ipopt





function create_tax_problem(n::Int64)
    A = 1:n;
    B = 1:n;
    T = [(i,j) for i in A for j in B];

    lambda = ones(n,n);
    w = wmin = 1;
    wmax = wmin + n - 1;
    w = [wmin + ((wmax-wmin)/(n-1))*(i-1) for i = 1:n]
    mu = [1.0 + i/n for i = 1:n]
    gamma = 1.0

    mu1 = mu + 1

    m = Model(solver = IpoptSolver())
    @variable(m, c[i in A,j in B] >= 0.0)# for (i,j) in T)
    @variable(m, y[i in A,j in B] >= 0.0)# for (i,j) in T )

    @NLobjective(m, Min, -sum( lambda[i,j] * (log(c[i,j]) - (y[i,j]/w[i])^mu[j] / mu[j])  - 0.5 * gamma * (c[i,j]^2 + y[i,j]^2) for (i,j) = T) )

    for (i,j) in T
      for (p,q) in T
        if i != p || j != p
          @NLconstraint(m,(log(c[i,j]) - (y[i,j]/w[i])^mu1[j] / mu[j]) - (log(c[p,q]) - (y[p,q]/w[i])^mu1[j] / mu[j]) >= 0)
        end
      end
    end

    @NLconstraint(m, sum(lambda[i] * (y[i,j] - c[i,j]) for (i,j) in T)  >= 0.0);

    return m
end
#println("Objective value: ", getobjectivevalue(m))
#println("x = ", getvalue(x))
#println("y = ", getvalue(y))



#status = solve(m)

#getvalue(x), getvalue(y)
#AbstractNLPModel


function one_phase_solve(m)
    nlp_raw = MathProgNLPModel(m);

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

    @assert(is_feasible(intial_it, my_par.comp_feas))
    iter, status, hist, t, err = one_phase_IPM(intial_it, my_par, timer);

    pause_advanced_timer(timer)

    print_timer_stats(timer)
end

srand(1)
m = create_tax_problem(5)
one_phase_solve(m)

#=srand(1)
m = create_tax_problem(5)
one_phase_solve(m)

srand(1)
m = create_tax_problem(10)
one_phase_solve(m)


srand(1)
m = create_tax_problem(20)
one_phase_solve(m)=#

ORG_STDOUT = STDOUT
file = open("test.txt", "w")
redirect_stdout(file)

srand(1)
m = create_tax_problem(10)
one_phase_solve(m)

STDOUT = ORG_STDOUT

close(file)


#srand(1)
#m = create_tax_problem(40)
#one_phase_solve(m)


#=srand(1)
m = create_tax_problem(5)
status = solve(m)

srand(1)
m = create_tax_problem(10)
status = solve(m)

srand(1)
m = create_tax_problem(20)
status = solve(m)=#

#srand(1)
#m = create_tax_problem(40)
#status = solve(m)

#finalize(nlp_raw)
