include("../include.jl")

using JuMP
using Ipopt


function formulation_3D(lambda,w,mu,mu1, A, B, C, alpha)
    T = [(i,j,k) for i in A for j in B for k in C];

    m = Model(solver = IpoptSolver())
    @variable(m, c[i in A,j in B,k in C] >= alpha[k])# for (i,j) in T)
    @variable(m, y[i in A,j in B,k in C] >= alpha[k])# for (i,j) in T )


    @NLobjective(m, Min, -sum( lambda[i,j,k] * (log(c[i,j,k]) - (y[i,j,k]/w[i])^mu1[j] / mu1[j]) for (i,j,k) = T) )

    for (i,j,k) in T
      for (p,q,r) in T
        if i != p || j != q || k != r
          @NLconstraint(m,(log(c[i,j,k]) - (y[i,j,k]/w[i])^mu1[j] / mu1[j]) - (log(c[p,q,r]) - (y[p,q,r]/w[i])^mu1[j] / mu1[j]) >= 0.0)
        end
      end
    end

    @NLconstraint(m, sum(lambda[i,j,k] * (y[i,j,k] - c[i,j,k]) for (i,j,k) in T)  >= 0.0);
    return m
end

na = 7
nb = 3
nc = 3

A = 1:na;
B = 1:nb;
C = 1:nc;

lambda = ones(na,nb,nc);
wmin = 2;
#wmax = wmin + na - 1;
wmax = 4
w = [wmin + ((wmax-wmin)/(na-1))*(i-1) for i = 1:na]
#mu = [1.0, 2.0, 3.0, 5.0, 8.0]
#mu = [1.0, 2.0, 4.0] #, 5.0, 8.0]

mu = [0.5, 1.0, 2.0]

gamma = 0.0

alpha = [0.0, 1.0, 1.5]

mu1 = mu + 1.0


srand(1)
m = formulation_3D(lambda,w,mu,mu1,A,B,C, alpha)
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
