mutable struct linear_solver_JULIA <: abstract_linear_system_solver
	_factor # TO DO, give type
	_factor_defined::Bool

	# options
	sym::Symbol # this determines if we use LU, cholesky or LDL.
	safe_mode::Bool
	recycle::Bool


  function linear_solver_JULIA(sym::Symbol, safe_mode::Bool, recycle::Bool)
			this = new();
			this.sym = sym
			this.safe_mode = safe_mode
			this.recycle = recycle
			this._factor_defined = false
			return this
  end
end

function ls_factor!(solver::linear_solver_JULIA, SparseMatrix::SparseMatrixCSC{Float64,Int64}, n::Int64, m::Int64, timer::class_advanced_timer)
			inertia_status_val = 1;

			start_advanced_timer(timer, "JULIA/factorize")
			if solver.sym == :unsymmetric
				# incomplete ---- how does one compute inertia????
				 solver._factor = lufact(SparseMatrix);
			elseif solver.sym == :definite
				# do a cholesky factorization
				@assert(m == 0)
				try
					if !solver.recycle || !solver._factor_defined
						#solver._factor = cholfact(Symmetric(SparseMatrix,:L));
						solver._factor = cholesky(Symmetric(SparseMatrix,:L));
						solver._factor_defined = true
					else
						cholfact!(solver._factor,Symmetric(SparseMatrix,:L));
					end
				catch(e)
					if typeof(e) == LinearAlgebra.PosDefException
						inertia_status_val = 0
					else
							println("ERROR in linear_solver_JULIA.ls_factor!")
							throw(e)
					end
				end
			elseif solver.sym == :symmetric
				# LDL factorization
				try
					if !solver.recycle || !solver._factor_defined
						#solver._factor = ldltfact(Symmetric(SparseMatrix,:L));
						solver._factor = ldlt(Symmetric(SparseMatrix,:L));
						solver._factor_defined = true
					else
						ldltfact!(solver._factor,Symmetric(SparseMatrix,:L));
					end
				catch(e)
					#PosDefException was added when updating Julia from 0.6.4 to 0.7.0 by Fadi Hamad
					#Because ArgumentError is not thrown anymore and it is being replaced by PosDefException
					#This is related to (PosDefException: matrix is not positive definite; Cholesky factorization failed.)
					if typeof(e) == ArgumentError || typeof(e) == LinearAlgebra.PosDefException || typeof(e) == ZeroPivotException
					#if typeof(e) == LinearAlgebra.PosDefException
							inertia_status_val = 0
					else
							println("ERROR in saddle_solver_JULIA.ls_factor!")
                                                        #println("###############################", typeof(e))
							throw(e)
					end
				end

				if inertia_status_val == 1
					# do something !!!!!!
					di = diag(solver._factor)
					tol = 1e-20
					pos_eigs =  sum(di .> tol)
					zero_eigs = sum((-tol .<= di) .& (di .<= tol))
					neg_eigs = sum(di .< -tol)
					nan_eigs = sum(isnan.(di))
					inf_eigs = sum(isinf.(di))
					if nan_eigs + inf_eigs <= 0
						inertia_status_val = inertia_status(pos_eigs, neg_eigs, zero_eigs, n, m)
					else
						if inf_eigs > 0
							warn("Inf appears in diagonal")
						end
						if nan_eigs > 0
							warn("NaN appears in diagonal")
						end
						inertia_status_val = 0
					end
				end
			else
				error("this.options.sym = " * string(solver.sym) * " not supported")
			end
			pause_advanced_timer(timer, "JULIA/factorize")

			return inertia_status_val; # inertia
end

function ls_solve!(solver::linear_solver_JULIA, my_rhs::Array{Float64,1}, my_sol::Array{Float64,1}, timer::class_advanced_timer)
	start_advanced_timer(timer, "JULIA/ls_solve")
  my_sol[1:length(my_sol)] = solver._factor \ my_rhs; #::UmfpackLU{Float64,Int64}
	pause_advanced_timer(timer, "JULIA/ls_solve")
end

function ls_solve(solver::linear_solver_JULIA, my_rhs::AbstractArray, timer::class_advanced_timer)
    start_advanced_timer(timer, "JULIA/ls_solve")
    #Vector was added when updating from Julia 0.6.4 to Julia 0.7.0 by Fadi Hamad
    #This is becasue my_rhs in somecases is being a SparseVector instead of being array.
    #Due to the code in shur.jl jac_res = eval_jac_T_prod and hess_res = hess_product are returning SparseMatrix now instead of Array
    sol = solver._factor \ Vector(my_rhs);
    pause_advanced_timer(timer, "JULIA/ls_solve")
    return sol
end
