include("../include.jl")

using JuMP
using Ipopt

function create_tax_problem(na::Int64)
    A = 1:na;

    mu = ones(na);
    lambda = ones(na)

    wmin = 1;
    wmax = wmin + na - 1;
    w = [wmin + ((wmax-wmin)/(na-1))*(i-1) for i = 1:na]


    m = Model(solver = IpoptSolver(print_level=5))
    @variable(m, c[i in A] >= 0.0)# for (i,j) in T)
    @variable(m, y[i in A] >= 0.0)# for (i,j) in T )

    @NLobjective(m, Min, -sum( lambda[i] * (log(c[i]) - (y[i]/w[i])^(mu[i]+1) / (mu[i]+1)) for i in A) )

    # Incentive
    for i in A
      for p in A
        if i != p
          @NLconstraint(m,(log(c[i]) - (y[i]/w[i])^(mu[i]+1) / (mu[i]+1))
        - (log(c[p]) - (y[p]/w[i])^(mu[i]+1) / (mu[i]+1)) >= 0.0)
        end
      end
    end

    @NLconstraint(m, sum(lambda[i] * (y[i] - c[i]) for i in A)  >= 0.0);

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

    return iter
end

m = create_tax_problem(40);
status = solve(m)
iter = one_phase_solve(m);

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
m = create_tax_problem(30)
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
