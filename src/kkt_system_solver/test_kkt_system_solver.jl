include("../linear_system_solvers/linear_system_solvers.jl")

include("kkt_system_solver.jl")

using advanced_timer

start_advanced_timer()
iter = toy_LP2()
dir_mu = 0.5

kkt_solver = Schur_KKT_solver()
kkt_solver.ls_solver = linear_solver_JULIA(:unsymmetric)
initialize!(kkt_solver)
factor!(kkt_solver, iter)
rhs = System_rhs()
rhs.primal_r = zeros(ncon(iter))
rhs.dual_r = zeros(dim(iter))
rhs.comp_r = zeros(ncon(iter))

compute_direction(kkt_solver, rhs, dir_mu)

kkt_solver = Schur_KKT_solver()
kkt_solver.ls_solver = linear_solver_MUMPS(:unsymmetric)
initialize!(kkt_solver)
factor!(kkt_solver, iter)
rhs = System_rhs()
rhs.primal_r = zeros(ncon(iter))
rhs.dual_r = zeros(dim(iter))
rhs.comp_r = zeros(ncon(iter))

compute_direction(kkt_solver, rhs, dir_mu)
