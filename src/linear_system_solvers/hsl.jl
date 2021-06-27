#println("HSL library not working")
print("loading HSL lib ... ")
using HSL
println("HSL lib loaded.")

mutable struct linear_solver_HSL <: abstract_linear_system_solver
	_factor # TO DO, give type
	#_factor_defined::Bool

	# options
	sym::Symbol # this determines if we use LU, cholesky or LDL.
	#safe_mode::Bool
	#recycle::Bool
	u::Float64


  function linear_solver_HSL(sym::Symbol, safe_mode::Bool, recycle::Bool, u::Float64)
			this = new();
			this.sym = sym
			this.u = u
			#this.safe_mode = safe_mode
			#this.recycle = recycle
			#this._factor_defined = false
			return this
  end
end

function ls_factor!(solver::linear_solver_HSL, SparseMatrix::SparseMatrixCSC{Float64,Int64}, n::Int64, m::Int64, timer::class_advanced_timer)
			inertia_status_val = 1;
            # 1 indicates inertia value is correct
            # any other value indicates it is incorrect

			start_advanced_timer(timer, "JULIA/factorize")
			if solver.sym == :symmetric
				#@assert(m == 0)
				solver._factor = Ma97(SparseMatrix,print_level=-1)
				solver._factor.control.u = solver.u
                ma97_factorize!(solver._factor);
                info = solver._factor.info
                if info.num_neg > m
                    inertia_status_val = 0
                elseif info.matrix_rank != n + m
                    inertia_status_val = 0
                end
			else
				error("this.options.sym = " * string(solver.sym) * " not supported")
			end
			pause_advanced_timer(timer, "JULIA/factorize")

			return inertia_status_val; # inertia
end

function ls_solve!(solver::linear_solver_HSL, my_rhs::Array{Float64,1}, my_sol::Array{Float64,1}, timer::class_advanced_timer)
	start_advanced_timer(timer, "JULIA/ls_solve")
    my_sol[1:length(my_sol)] = ma97_solve(solver._factor, my_rhs);
	pause_advanced_timer(timer, "JULIA/ls_solve")
end

function ls_solve(solver::linear_solver_HSL, my_rhs::AbstractArray, timer::class_advanced_timer)
	start_advanced_timer(timer, "JULIA/ls_solve")
    #println("---------------------------------------typeof(my_rhs): ", typeof(my_rhs))
    #println("-----------------------------------------------my_rhs: ", my_rhs)
    if typeof(my_rhs) == SparseVector{Float64, Int64}
         my_rhs = Vector(my_rhs)
    end
    sol = ma97_solve(solver._factor, my_rhs);
	pause_advanced_timer(timer, "JULIA/ls_solve")

	return sol
end

