This page describes each of the columns in the solver output.

it = iteration number
step = step type, either stabilization step (s) or aggressive step (a)
eta = *targeted* reduction in feasibility/barrier parameter, eta=1 for stabilization steps eta<1 for aggressive steps
α_P = primal step size
α_D = dual step size
ls = number of points trialled during line search
|dx| = infinity norm of primal direction size
|dy| = infinity norm of dual direction size
N err = relative error in linear system solves.
mu = value of barrier parameter
dual = gradient of lagragian scaled by largest dual variable
primal = error in primal feasibility
cmp scaled = \| Sy \|/(1 + \|y\|)
inf = how close to infeasible problem is, values close to zero indicate infeasibility
delta = size of peturbuation
\#fac  = number of factorizations (split into two numbers -- first is how many factorization needed to ensure primal schur complement is positive definite, second number represents total number of factorizations including any increases in delta to avoid issues when direction quality is very poor)
|x| = infinity norm of x variables
|y| = infinity norm of y variables
∇phi = gradient of log barrier
phi = value of log barrier
