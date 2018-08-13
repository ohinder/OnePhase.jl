println("Loading linear_system_solvers ... ")
@compat abstract type abstract_linear_system_solver end

## define linear system solvers.
## note that linear system solvers need to be coupled with a KKT system solver
## these solvers can be found in the `kkt_system_solver` folder and
## reduces the linear system that needs to be solved, e.g., does a "primal schur complement"

include("julia.jl")
USE_HSL = true
if USE_HSL
	try
		include("hsl.jl")
	catch (e)
		warn(e)
	end
end
#include("matlab.jl")
if USE_MUMPS
	include("mumps_wrapper.jl")
end

function initialize!(solver::abstract_linear_system_solver)

end

function finalize!(solver::abstract_linear_system_solver)

end

function inertia_status(pos_eigs::Int64, neg_eigs::Int64, zero_eigs::Int64, num_vars::Int64, num_constraints::Int64)
	try
		#println("inertia_status called") # ????
		@assert(pos_eigs + neg_eigs + zero_eigs == num_vars + num_constraints)
		# return number instead of true or false???

		if pos_eigs == num_vars && neg_eigs == num_constraints
			return true
		elseif pos_eigs > num_vars || neg_eigs + pos_eigs != num_vars + num_constraints
			# Warning: numerical instability in LDL factorization
			if pos_eigs > num_vars
				# more positive eigenvalues than variables
				return false
			else
				# zero eigenvalues
				return false
			end
		else
			return false
		end
	catch e
		println("ERROR in compute_inertia")
		throw(e)
	end
end
