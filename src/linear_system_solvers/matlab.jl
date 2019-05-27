using MATLAB

type linear_solver_MATLAB <: abstract_linear_system_solver
  _SparseMatrix::SparseMatrixCSC{Float64,Int64}

	# options
	sym::Symbol

  function linear_solver_MATLAB()
      # start matlab session
      @matlab begin
      end

      return new();
  end
end

function ls_factor!(solver::linear_solver_MATLAB, n::Int64, m::Int64)
    # solver = linear system solver which calls matlab to solve linear systems
    # n = number of variables
    # m = number of constraints

			sparse_matrix = solver. _SparseMatrix
			@mput sparse_matrix

			@matlab begin
				mat_L, mat_D, mat_P, mat_S = ldl(sparse_matrix,0.5);
			end

			@mget mat_D

			pos_eigs, neg_eigs, zero_eigs = eigsigns_of_B(mat_D, 0.0);
			return inertia_status(pos_eigs, neg_eigs, zero_eigs, n, m)
end

function ls_solve(solver::linear_solver_MATLAB, my_rhs::AbstractArray)
      @mput my_rhs

      @matlab begin
        mat_solution = mat_S * ( mat_P* (mat_L'\ (mat_D\ (mat_L\ (mat_P' * ( mat_S * my_rhs ) ) ) ) ) );
      end

      @mget mat_solution

      return mat_solution;
end
