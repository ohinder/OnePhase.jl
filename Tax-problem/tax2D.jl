include("../include.jl")

using JuMP
using Ipopt


function formulation1(lambda,w,mu,mu1,A,B,T, lb_var)
    m = Model(solver = IpoptSolver())
    @variable(m, c[i in A,j in B] >= lb_var)# for (i,j) in T)
    @variable(m, y[i in A,j in B] >= lb_var)# for (i,j) in T )


    @NLobjective(m, Min, -sum( lambda[i,j] * (log(c[i,j]) - (y[i,j]/w[i])^mu1[j] / mu1[j]) for (i,j) = T) )

    for (i,j) in T
      for (p,q) in T
        if i != p || j != p
          @NLconstraint(m,(log(c[i,j]) - (y[i,j]/w[i])^mu1[j] / mu1[j]) - (log(c[p,q]) - (y[p,q]/w[i])^mu1[j] / mu1[j]) >= 0.0)
        end
      end
    end

    @NLconstraint(m, sum(lambda[i] * (y[i,j] - c[i,j]) for (i,j) in T)  >= 0.0);
    return m
end

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
    #intial_it.point.mu *= 10.0
    #intial_it.point.y *= 10.0
    pause_advanced_timer(timer, "INIT")

    pause_advanced_timer(timer)
    print_timer_stats(timer)

    start_advanced_timer(timer)

    @assert(is_feasible(intial_it, my_par.comp_feas))
    iter, status, hist, t, err = one_phase_IPM(intial_it, my_par, timer);

    pause_advanced_timer(timer)

    print_timer_stats(timer)

    return iter
end

na = 17
nb = 5

A = 1:na;
B = 1:nb;
T = [(i,j) for i in A for j in B];

lambda = ones(na,nb);
wmin = 1;
#wmax = wmin + na - 1;
wmax = 5
w = [wmin + ((wmax-wmin)/(na-1))*(i-1) for i = 1:na]
#mu = [1.0, 2.0, 3.0, 5.0, 8.0]
mu = [1.0, 1.5, 2.0, 2.5, 3.0]

gamma = 0.0

mu1 = mu + 1.0

#lb = 0.0
lb_var = 0.1

srand(1)
m = formulation1(lambda,w,mu,mu1,A,B,T, lb_var)
iter = one_phase_solve(m);

#=
ORG_STDOUT = STDOUT
file = open("test.txt", "w")
redirect_stdout(file)

srand(1)
m = create_tax_problem(10)
one_phase_solve(m)

STDOUT = ORG_STDOUT

close(file)=#
