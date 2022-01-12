using LinearAlgebra

################################################################################
## Defines linear system solvers.
## Note that linear system solvers need to be coupled with a KKT system solver
## these solvers can be found in the `kkt_system_solver` folder and
## reduces the linear system that needs to be solved, e.g., does a "primal schur complement"
################################################################################

println("Loading linear_system_solvers ... ")
abstract type abstract_linear_system_solver end

include("julia.jl")
global USE_HSL = true
setUSE_HSL(use_hsl) = (global USE_HSL = use_hsl)

export USE_HSL, setUSE_HSL, loadHSL

function loadHSL(hsl_dir)
    if USE_HSL
	try
		hsl_code_location = string(hsl_dir, "hsl.jl")
		include(hsl_code_location)
	catch (e)
		println("Loading HSL failed:")
		@warn(e)
		println("Continuing although you will not be able to choose HSL as a linear solver ...")
	end
    end
    
end

loadHSL("./")

#include("matlab.jl")
if USE_MUMPS
	include("mumps_wrapper.jl")
end

function initialize!(solver::abstract_linear_system_solver)

end

function finalize!(solver::abstract_linear_system_solver)

end

function inertia_status(pos_eigs::Int64, neg_eigs::Int64, zero_eigs::Int64, num_vars::Int64, num_constraints::Int64)
	###########################################################################################
	# INPUT:
	# pos_eigs, neg_eigs, zero_eigs = number of positive, negative and zero eigenvalues in the matrix:
	# [[ H A'];
	# [ A D ]];
	# num_vars = number of variables in problem
	# num_constraints = number of constraints in problem
	# OUTPUT:
	# is the inertia good, i.e., is H + A' D^{-1} A positive semidefinite?
	# this ensures that our direction will be a descent direction on the shifted log barrier.
	###########################################################################################

	try
		#println("inertia_status called") # ????
		if(pos_eigs + neg_eigs + zero_eigs != num_vars + num_constraints)
			println("pos_eigs = $pos_eigs")
			println("neg_eigs = $neg_eigs")
			println("zero_eigs = $pos_eigs")
			println("num_vars = $num_vars")
			println("num_constraints = $num_constraints")
			error("pos_eigs + neg_eigs + zero_eigs != num_vars + num_constraints")
		end
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
