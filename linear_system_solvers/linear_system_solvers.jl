println("Loading linear_system_solvers ...")
abstract abstract_linear_system_solver;

include("julia.jl")
#include("matlab.jl")
if MUMPS
	include("mumps_wrapper.jl")
end

function initialize!(solver::abstract_linear_system_solver)

end

function finalize!(solver::abstract_linear_system_solver)

end

function inertia_status(pos_eigs::Int64, neg_eigs::Int64, zero_eigs::Int64, num_vars::Int64, num_constraints::Int64)
	try
		@assert(pos_eigs + neg_eigs + zero_eigs == num_vars + num_constraints)
		#println(pos_eigs, " ", neg_eigs);

		if pos_eigs == num_vars && neg_eigs == num_constraints
			return true
		elseif pos_eigs > num_vars || neg_eigs + pos_eigs != num_vars + num_constraints
			println("Warning: numerical instability in LDL factorization")
			if pos_eigs > num_vars
				println("more positive eigenvalues than variables")
				return false
			else
				println("zero eigenvalues")
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
