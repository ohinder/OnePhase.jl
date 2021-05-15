include("linear_system_solvers.jl")
using advanced_timer
#using SparseArrays

start_advanced_timer()
solver_julia = linear_solver_JULIA(:unsymmetric)
initialize!(solver_julia)
ls_factor!(solver_julia, SparseArrays.speye(10), 10, 0)
ls_solve(solver_julia, rand(10))
#mumps_sym(:symmetric)

solver_mumps = linear_solver_MUMPS(:unsymmetric)
initialize!(solver_mumps)
ls_factor!(solver_mumps, SparseArrays.speye(10), 10, 0)
ls_solve(solver_mumps, rand(10))
pause_advanced_timer()



## attempt to produce minimal example that re-produces seg-faults
##

using MUMPS, MPI
MPI.Init()

cntl = MUMPS.default_cntl64[:]; # copy
cntl[1] = 0.01

icntl = MUMPS.default_icntl[:]; # copy

icntl[4] = 1;
icntl[10] = 2; # 2 iterative refinement steps
icntl[14] = 200.0 #1000.0;

n = 1000

my_kkt_solver = Symmetric_KKT_solver()
my_kkt_solver.ls_solver = linear_solver_MUMPS(:symmetric)

my_kkt_solver.ls_solver._factor = MUMPS.Mumps{Float64}(MUMPS.mumps_symmetric, icntl, cntl)

function solves!(my_kkt_solver::abstract_KKT_system_solver)
  MUMPS.associate_rhs!(my_kkt_solver.ls_solver._factor, rand(n));
  MUMPS.solve!(my_kkt_solver.ls_solver._factor);
  sol = MUMPS.get_solution(my_kkt_solver.ls_solver._factor);
  @show LinearAlgebra.norm(sol)
end

function fact!(my_kkt_solver::abstract_KKT_system_solver)
    MUMPS.associate_matrix!(my_kkt_solver.ls_solver._factor, SparseArrays.sparse(rand(n,n)));
    MUMPS.factorize!(my_kkt_solver.ls_solver._factor);
end

for i = 1:50
    fact!(my_kkt_solver)
    for i = 1:50
      solves!(my_kkt_solver)
    end
end
