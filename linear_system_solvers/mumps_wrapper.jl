println("loading mumps lib ...")
using MPI
using MUMPS # needs to be intialized before Ipopt don't understand why
println("mumps lib loaded.")

type linear_solver_MUMPS <: abstract_linear_system_solver
	_factor::MUMPS.Mumps{Float64}
	M::SparseMatrixCSC{Float64,Int64}
	sym::Symbol
	safe_mode::Bool

  function linear_solver_MUMPS(sym::Symbol, safe_mode::Bool)
      this = new();
			this.sym = sym
			this.safe_mode = safe_mode

			return this
  end
end

# intialize and finalize function???
# MPI ???
function mumps_sym(sym::Symbol)
	if sym == :symmetric
			return MUMPS.mumps_symmetric
	elseif sym == :unsymmetric
			return MUMPS.mumps_unsymmetric
	elseif sym == :definite
			return MUMPS.mumps_definite
	else
			error("this symmetry symbol is not understood")
	end
end

function create_mumps_factor(solver::linear_solver_MUMPS)
	cntl = MUMPS.default_cntl64[:]; # copy
	cntl[1] = 0.01

	icntl = MUMPS.default_icntl[:]; # copy

	icntl[4] = 1;
	icntl[10] = 1; # 2 iterative refinement steps
	icntl[14] = 200.0 #1000.0;


	return MUMPS.Mumps{Float64}(mumps_sym(solver.sym), icntl, cntl	);  # Real, general unsymmetric
end

function create_mumps_factor2(solver::linear_solver_MUMPS)
	icntl = get_icntl(verbose=false);
	return MUMPS.Mumps{Float64}(mumps_sym(solver.sym), icntl, MUMPS.default_cntl64);  # Real, general unsymmetric
end

function initialize!(solver::linear_solver_MUMPS)
		if ~MPI.Initialized()
			MPI.Init()
		end

    solver._factor = create_mumps_factor2(solver)
end

function finalize!(solver::linear_solver_MUMPS)
    MUMPS.finalize(solver._factor);
end

function ls_factor!(solver::linear_solver_MUMPS, SparseMatrix::SparseMatrixCSC{Float64,Int64}, n::Int64, m::Int64, timer::class_advanced_timer)
		#if
		return ls_factor_LDL!(solver, SparseMatrix, n, m, timer)
end

# A-matrix
function ls_factor_LDL!(solver::linear_solver_MUMPS, SparseMatrix::SparseMatrixCSC{Float64,Int64}, n::Int64, m::Int64, timer::class_advanced_timer)
		start_advanced_timer(timer,"MUMPS")
		@assert(size(SparseMatrix,1) == n + m)
		@assert(size(SparseMatrix,2) == n + m)

		solver.M = deepcopy(SparseMatrix)

		start_advanced_timer(timer,"MUMPS/associate_matrix")
    associate_matrix!(solver._factor, solver.M);
		pause_advanced_timer(timer,"MUMPS/associate_matrix")

		start_advanced_timer(timer,"MUMPS/factorize")
		factorize!(solver._factor);
		pause_advanced_timer(timer,"MUMPS/factorize")

		zero_eigs = solver._factor.infog[1] == -10 || solver._factor.infog[1] == -2

		pause_advanced_timer(timer,"MUMPS")

		if solver._factor.infog[12] > m
				return -1
		end

		if zero_eigs
				return 0
		end

		return 1
end

function ls_solve(solver::linear_solver_MUMPS, my_rhs::AbstractArray, timer::class_advanced_timer)
		start_advanced_timer(timer,"MUMPS")
		start_advanced_timer(timer,"MUMPS/solve")

		if solver.safe_mode
			MUMPS.associate_matrix!(solver._factor, solver.M);
			MUMPS.factorize!(solver._factor);
		end

		associate_rhs!(solver._factor, reshape(my_rhs,length(my_rhs),1))
		solve!(solver._factor)
		sol = deepcopy(get_solution(solver._factor))
		#sol = solve(solver._factor, my_rhs)

		pause_advanced_timer(timer,"MUMPS/solve")
		pause_advanced_timer(timer,"MUMPS")

    return sol[:]
end
