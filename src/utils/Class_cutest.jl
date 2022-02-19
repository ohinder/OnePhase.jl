#importall CUTEst
#using SparseArrays
#import CUTEst
include("../parameters.jl")
include("../linear_system_solvers/linear_system_solvers.jl")
include("../kkt_system_solver/include.jl")
include("../line_search/line_search.jl")
include("../IPM/ipm.jl")
include("../JuMPinterface.jl")
export get_original_x, get_y, number_constraints

mutable struct Class_bounds
    l_i::Array{Int64,1}
    u_i::Array{Int64,1}
    l::Array{Float64,1}
    u::Array{Float64,1}

    function Class_bounds(lb_vec::Array{Float64,1}, ub_vec::Array{Float64,1})
        @assert(length(lb_vec) == length(ub_vec))
        this = new([],[],[],[]);
        for i = 1:length(lb_vec)
            if lb_vec[i] > -Inf
               push!(this.l_i, i)
               push!(this.l, lb_vec[i])
            end

            if ub_vec[i] < Inf
               push!(this.u_i, i)
               push!(this.u, ub_vec[i])
            end
        end

        return this;
    end
end

function _i_not_fixed(m::NLPModels.AbstractNLPModel)
    return (1:m.meta.nvar)[m.meta.lvar .!= m.meta.uvar]
end

function _i_not_fixed(variable_info::Vector{VariableInfo})
    vector_indixes = zeros(Int64, 0)
    for i in 1:length(variable_info)
        if !variable_info[i].is_fixed
            push!(vector_indixes, i)
        end
    end
    return vector_indixes
end

function extract_lcon_ucon(solver::OnePhaseSolver)
    lcon = zeros(Float64, 0)
    ucon = zeros(Float64, 0)

    #Linear constraints
    for i in 1:length(solver.linear_eq_constraints)
        lower = solver.linear_eq_constraints[i].set.value
        push!(lcon, lower)
        push!(ucon, lower)
    end

    for i in 1:length(solver.linear_ge_constraints)
        lower = solver.linear_ge_constraints[i].set.lower
        upper = Inf
        push!(lcon, lower)
        push!(ucon, upper)
    end

    for i in 1:length(solver.linear_le_constraints)
        lower = -Inf
        upper = solver.linear_le_constraints[i].set.upper
        push!(lcon, lower)
        push!(ucon, upper)
    end

    for i in 1:length(solver.linear_int_constraints)
        lower = solver.linear_int_constraints[i].set.lower
        upper = solver.linear_int_constraints[i].set.upper
        push!(lcon, lower)
        push!(ucon, upper)
    end

    #Quadratic Constraints
    for i in 1:length(solver.quadratic_eq_constraints)
        lower = solver.quadratic_eq_constraints[i].set.lower
        push!(lcon, lower)
        push!(ucon, lower)
    end

    for i in 1:length(solver.quadratic_ge_constraints)
        lower = solver.quadratic_ge_constraints[i].set.lower
        upper = Inf
        push!(lcon, lower)
        push!(ucon, upper)
    end

    for i in 1:length(solver.quadratic_le_constraints)
        lower = -Inf
        upper = solver.quadratic_le_constraints[i].set.upper
        push!(lcon, lower)
        push!(ucon, upper)
    end

    for i in 1:length(solver.quadratic_int_constraints)
        lower = solver.quadratic_int_constraints[i].set.lower
        upper = solver.quadratic_int_constraints[i].set.upper
        push!(lcon, lower)
        push!(ucon, upper)
    end

    #Nonlinear Constraints
    for i in 1:length(solver.nlp_data.constraint_bounds)
        push!(lcon, solver.nlp_data.constraint_bounds[i].lower)
        push!(ucon, solver.nlp_data.constraint_bounds[i].upper)
    end

    return lcon, ucon
end

function extract_lvar_uvar(variable_info::Vector{VariableInfo}, ind::Vector{Int64})
    lvar = zeros(Float64, 0)
    uvar = zeros(Float64, 0)
    for i in ind
        push!(lvar, variable_info[i].lower_bound)
        push!(uvar, variable_info[i].upper_bound)
    end
    return lvar, uvar
end

mutable struct Class_CUTEst <: abstract_nlp
    nlp::Union{NLPModels.AbstractNLPModel, Nothing}
    solver::Union{OnePhaseSolver, Nothing}

    bcon::Class_bounds
    bvar::Class_bounds

    function Class_CUTEst(nlp::NLPModels.AbstractNLPModel)
        ind = _i_not_fixed(nlp)
        return new(nlp, nothing, Class_bounds(nlp.meta.lcon[:], nlp.meta.ucon[:]), Class_bounds(nlp.meta.lvar[ind], nlp.meta.uvar[ind]))
    end

    function Class_CUTEst(solver::OnePhaseSolver)
        ind = _i_not_fixed(solver.variable_info)
        lcon, ucon = extract_lcon_ucon(solver)
        lvar, uvar = extract_lvar_uvar(solver.variable_info, ind)
        return new(nothing, solver, Class_bounds(lcon, ucon), Class_bounds(lvar, uvar))
    end
end

function suggested_starting_point(m::Class_CUTEst)
    if m.nlp != nothing
        ind = _i_not_fixed(m.nlp)
        return deepcopy(m.nlp.meta.x0[ind])
    end

    ind = _i_not_fixed(m.solver.variable_info)
    non_fixed_variables = m.solver.variable_info[ind]
    starting_points_vector = zeros(0)
    for variable in non_fixed_variables
        push!(starting_points_vector, variable.start == nothing ? 0.0 : variable.start)
    end

    return deepcopy(starting_points_vector)
end

#IMPORTANT_
# function suggested_starting_point(m::Class_CUTEst)
#     if m.nlp != nothing
#         ind = _i_not_fixed(m.nlp)
#         return deepcopy(m.nlp.meta.x0[ind])
#     end
#     ind = _i_not_fixed(m.solver.variable_info)
#     non_fixed_variables = m.solver.variable_info[ind]
#     starting_points_vector = zeros(0)
#     for variable in non_fixed_variables
#         push!(starting_points_vector, variable.start == nothing ? 0.0 : variable.start)
#     end
#     return deepcopy(starting_points_vector)
# end

function ncon(m::Class_CUTEst)
    return nbounds_orginal(m) + ncons_orginal(m)
end

function lb(x::Array{Float64,1}, bd::Class_bounds)
    return x[bd.l_i] - bd.l
end

function ub(x::Array{Float64,1}, bd::Class_bounds)
    return bd.u - x[bd.u_i]
end

function number_constraints(solver::OnePhaseSolver)
    linear_constraints_count = number_constraints_linear(solver)
    quadratic_constraints_count = number_constraints_quadratic(solver)
    non_linear_constraints_count = number_constraints_non_linear(solver)
    return linear_constraints_count + quadratic_constraints_count + non_linear_constraints_count
end

function number_constraints_linear(solver::OnePhaseSolver)
    less_than_constraints_count = length(solver.linear_le_constraints)
    greater_than_constraints_count = length(solver.linear_ge_constraints)
    equality_constraints_count = length(solver.linear_eq_constraints)
    interval_constraints_count = length(solver.linear_int_constraints)
    return less_than_constraints_count + greater_than_constraints_count + equality_constraints_count + interval_constraints_count
end

function number_constraints_quadratic(solver::OnePhaseSolver)
    less_than_constraints_count = length(solver.quadratic_le_constraints)
    greater_than_constraints_count = length(solver.quadratic_ge_constraints)
    equality_constraints_count = length(solver.quadratic_eq_constraints)
    interval_constraints_count = length(solver.quadratic_int_constraints)
    return less_than_constraints_count + greater_than_constraints_count + equality_constraints_count + interval_constraints_count
end

function number_constraints_non_linear(solver::OnePhaseSolver)
    return length(solver.nlp_data.constraint_bounds)
end

function linear_cons(m::Class_CUTEst)
    if m.nlp != nothing
        is_linear = zeros(m.nlp.meta.ncon)
        is_linear[m.nlp.meta.lin] .= 1.0
        vec = [is_linear[m.bcon.l_i]; is_linear[m.bcon.u_i]; ones(nbounds_orginal(m))]
        return 1.0 .== vec
    end

    is_linear = zeros(number_constraints(m.solver))
    is_linear[[i for i in 1:number_constraints_linear(m.solver)]] .= 1.0
    vec = [is_linear[m.bcon.l_i]; is_linear[m.bcon.u_i]; ones(nbounds_orginal(m))]
    return 1.0 .== vec
end

#IMPORTANT_
# function linear_cons(m::Class_CUTEst)
#     if m.nlp != nothing
#         is_linear = zeros(m.nlp.meta.ncon)
#         is_linear[m.nlp.meta.lin] .= 1.0
#         vec = [is_linear[m.bcon.l_i]; is_linear[m.bcon.u_i]; ones(nbounds_orginal(m))]
#         return 1.0 .== vec
#     end
#
#     is_linear = zeros(number_constraints(m.solver))
#
#     is_linear[[i for i in 1:number_constraints_linear(m.solver)]] .= 1.0
#     vec = [is_linear[m.bcon.l_i]; is_linear[m.bcon.u_i]; ones(nbounds_orginal(m))]
#
#     return 1.0 .== vec
# end

function ineq_cons(m::Class_CUTEst)
    if m.nlp != nothing
        is_ineq = zeros(m.nlp.meta.ncon)
        is_ineq[m.nlp.meta.lcon .== m.nlp.meta.ucon] .= 1.0
        vec = [is_ineq[m.bcon.l_i]; is_ineq[m.bcon.u_i]; ones(nbounds_orginal(m))]
        return 1.0 .== vec
    end

    solver = m.solver
    is_ineq =  zeros(number_constraints(solver))

    index = 1 + length(solver.linear_le_constraints) + length(solver.linear_ge_constraints)

    is_linear[[i for i in index:length(solver.linear_eq_constraints)]] .= 1.0

    index += length(solver.linear_eq_constraints)

    for i in 1:length(solver.linear_int_constraints)
        lower = solver.linear_int_constraints[i].set.lower
        upper = solver.linear_int_constraints[i].set.upper
        if lower == upper
            is_ineq[index] .= 1.0
        end
        index +=1
    end

    index += length(solver.quadratic_le_constraints) + length(solver.quadratic_ge_constraints)

    is_linear[[i for i in index:length(solver.quadratic_eq_constraints)]] .= 1.0

    index += length(solver.quadratic_eq_constraints)

    for i in 1:length(solver.quadratic_int_constraints)
        lower = solver.quadratic_int_constraints[i].set.lower
        upper = solver.quadratic_int_constraints[i].set.upper
        if lower .== upper
            is_ineq[index] .= 1.0
        end
        index +=1
    end

    #Nonlinear Constraints
    for i in 1:length(solver.nlp_data.constraint_bounds)
        lower = solver.nlp_data.constraint_bounds[i].lower
        upper = solver.nlp_data.constraint_bounds[i].upper
        if lower .== upper
            is_ineq[index] .= 1.0
        end
        index +=1
    end

    vec = [is_ineq[m.bcon.l_i]; is_ineq[m.bcon.u_i]; ones(nbounds_orginal(m))]
    return 1.0 .== vec
end

#IMPORTANT_
# function ineq_cons(m::Class_CUTEst)
#     if m.nlp != nothing
#         is_ineq = zeros(m.nlp.meta.ncon)
#         is_ineq[m.nlp.meta.lcon .== m.nlp.meta.ucon] .= 1.0
#         vec = [is_ineq[m.bcon.l_i]; is_ineq[m.bcon.u_i]; ones(nbounds_orginal(m))]
#
#         return 1.0 .== vec
#     end
#
#     solver = m.solver
#     is_ineq =  zeros(number_constraints(solver))
#
#     index = 1 + length(solver.linear_le_constraints) + length(solver.linear_ge_constraints)
#
#     is_linear[[i for i in index:length(solver.linear_eq_constraints)]] .= 1.0
#
#     index += length(solver.linear_eq_constraints)
#
#     for i in 1:length(solver.linear_int_constraints)
#         lower = solver.linear_int_constraints[i].set.lower
#         upper = solver.linear_int_constraints[i].set.upper
#         if lower == upper
#             is_ineq[index] .= 1.0
#         end
#         index +=1
#     end
#
#     index += length(solver.quadratic_le_constraints) + length(solver.quadratic_ge_constraints)
#
#     is_linear[[i for i in index:length(solver.quadratic_eq_constraints)]] .= 1.0
#
#     index += length(solver.quadratic_eq_constraints)
#
#     for i in 1:length(solver.quadratic_int_constraints)
#         lower = solver.quadratic_int_constraints[i].set.lower
#         upper = solver.quadratic_int_constraints[i].set.upper
#         if lower .== upper
#             is_ineq[index] .= 1.0
#         end
#         index +=1
#     end
#
#     #Nonlinear Constraints
#     for i in 1:length(solver.nlp_data.constraint_bounds)
#         lower = solver.nlp_data.constraint_bounds[i].lower
#         upper = solver.nlp_data.constraint_bounds[i].upper
#         if lower .== upper
#             is_ineq[index] .= 1.0
#         end
#         index +=1
#     end
#
#     vec = [is_ineq[m.bcon.l_i]; is_ineq[m.bcon.u_i]; ones(nbounds_orginal(m))]
#
#     return 1.0 .== vec
# end

function nbounds_orginal(nlp::Class_CUTEst)
    return length(nlp.bvar.l_i) + length(nlp.bvar.u_i)
end

function ncons_orginal(nlp::Class_CUTEst)
    return length(nlp.bcon.l_i) + length(nlp.bcon.u_i)
end

function cons_indicies(nlp::Class_CUTEst)
    m = ncons_orginal(nlp)
    if m > 0
      return 1:m
    else
      return []
    end
end

#function linear_indicies(nlp::Class_CUTEst)
#    #nlp.nlp.meta.lin
#end


function bound_indicies(nlp::Class_CUTEst)
    m = ncons_orginal(nlp)
    r = nbounds_orginal(nlp)
    if r > 0
      return (m + 1):(m + r)
    else
      return [];
    end
end

function y_l_con(y::Array{Float64,1}, m::Class_CUTEst)
    return y[1:length(m.bcon.l)]
end

function y_u_con(y::Array{Float64,1}, m::Class_CUTEst)
    n_lcon = length(m.bcon.l)
    return y[(n_lcon + 1):(n_lcon + length(m.bcon.u))]
end

function eval_f(m::Class_CUTEst, x::Array{Float64,1})
    if m.nlp != nothing
        return obj(m.nlp, _cute_x(m, x) )
    end

    return MOI.eval_objective(m.solver.nlp_data.evaluator, _cute_x(m, x))
end

#IMPORTANT_
# function eval_f(m::Class_CUTEst, x::Array{Float64,1})
#     if m.nlp != nothing
#         return obj(m.nlp, _cute_x(m, x) )
#     end
#
#     return MOI.eval_objective(m.solver.nlp_data.evaluator, _cute_x(m, x))
# end

function evaluateScalarAffineFunctionConstraints(constraints::Vector{Main.OnePhase.ConstraintInfo{MathOptInterface.ScalarAffineFunction{Float64}}}, x::Array{Float64,1})
    constraint_values = zeros(Float64, 0)
    for constraint in constraints
        value = constraint.func.constant
        for term in constraint.func.terms
            value += term.coefficient * x[term.variable.value]
        end
        push!(constraint_values, value)
    end
    return constraint_values
end

function evaluateScalarQuadraticFunctionConstraints(constraints::Vector{Main.OnePhase.ConstraintInfo{MathOptInterface.ScalarQuadraticFunction{Float64}}}, x::Array{Float64,1})
    constraint_values = zeros(Float64, 0)
    for constraint in constraints
        value = constraint.func.constant
        for term in constraint.func.affine_terms
            value += term.coefficient * x[term.variable.value]
        end
        for term in constraint.func.quadratic_terms
            value += term.coefficient * x[term.variable_1.value] * x[term.variable_2.value]
        end
        push!(constraint_values, value)
    end
    return constraint_values
end

function evaluateScalarAffineFunctionConstraintsJacobianStructure(constraints::Vector{Main.OnePhase.ConstraintInfo{MathOptInterface.ScalarAffineFunction{Float64}}}, x::Array{Float64,1})
    jac_structures = nothing
    count = 1
    for constraint in constraints
        jac_structure = zeros(Float64, length(x))
        for term in constraint.func.terms
            index = term.variable.value
            value = term.coefficient
            jac_structure[index] = value
        end
        if count == 1
            jac_structures = jac_structure
        else
            jac_structures = hcat(jac_structures, jac_structure)
        end
        count += 1
    end
    return jac_structures
end

function evaluateScalarQuadraticFunctionConstraintsJacobianStructure(constraints::Vector{Main.OnePhase.ConstraintInfo{MathOptInterface.ScalarQuadraticFunction{Float64}}}, x::Array{Float64,1})
    jac_structures = nothing
    count = 1
    for constraint in constraints
        jac_structure = zeros(Float64, length(x))
        for term in constraint.func.affine_terms
            index = term.variable.value
            value = term.coefficient
            jac_structure[index] = value
        end
        if count == 1
            jac_structures = jac_structure
        else
            jac_structures = hcat(jac_structures, jac_structure)
        end
        count += 1
    end
    return jac_structures
end

function eval_a(m::Class_CUTEst, x::Array{Float64,1})
    if m.nlp != nothing
        a = cons(m.nlp, _cute_x(m, x) )
        if a == Float64[]
          a = zeros(1)
        end
        a_1 = a
        return [lb(a, m.bcon); ub(a, m.bcon); lb(x, m.bvar); ub(x, m.bvar)];
    end

    solver = m.solver
    # linear_constraints = vcat(solver.linear_le_constraints, solver.linear_ge_constraints, solver.linear_eq_constraints, solver.linear_int_constraints)
    # linear_constraints_evaluation = evaluateScalarAffineFunctionConstraints(linear_constraints, x)
    #
    # quadratic_constraints = vcat(solver.quadratic_le_constraints, solver.quadratic_ge_constraints, solver.quadratic_eq_constraints, solver.quadratic_int_constraints)
    # quadratic_constraints_evaluation = evaluateScalarQuadraticFunctionConstraints(quadratic_constraints, x)

    linear_constraints = vcat(solver.linear_eq_constraints, solver.linear_ge_constraints, solver.linear_le_constraints, solver.linear_int_constraints)
    linear_constraints_evaluation = evaluateScalarAffineFunctionConstraints(linear_constraints, x)

    quadratic_constraints = vcat(solver.quadratic_eq_constraints, solver.quadratic_ge_constraints, solver.quadratic_le_constraints, solver.quadratic_int_constraints)
    quadratic_constraints_evaluation = evaluateScalarQuadraticFunctionConstraints(quadratic_constraints, x)

    non_linear_constraints_evaluation = Float64[]
    if typeof(m.solver.nlp_data.evaluator) != EmptyNLPEvaluator
        non_linear_constraints_evaluation = zeros(Float64, length(m.solver.nlp_data.evaluator.constraints))
        MOI.eval_constraint(m.solver.nlp_data.evaluator, non_linear_constraints_evaluation, x)
    end

    a = vcat(linear_constraints_evaluation, quadratic_constraints_evaluation, non_linear_constraints_evaluation)
    if a == Float64[]
      a = zeros(1)
    end

    return [lb(a, m.bcon); ub(a, m.bcon); lb(x, m.bvar); ub(x, m.bvar)];
end

#IMPORTANT_
# function eval_a(m::Class_CUTEst, x::Array{Float64,1})
#     if m.nlp != nothing
#         a = cons(m.nlp, _cute_x(m, x) )
#         if a == Float64[]
#           a = zeros(1)
#         end
#         return [lb(a, m.bcon); ub(a, m.bcon); lb(x, m.bvar); ub(x, m.bvar)];
#     end
#
#     # println("m.solver.linear_ge_constraints: ", m.solver.linear_ge_constraints[1].func)
#     # println("m.solver.linear_ge_constraints: ", m.solver.linear_ge_constraints[1].set)
#     solver = m.solver
#     linear_constraints = vcat(solver.linear_le_constraints, solver.linear_ge_constraints, solver.linear_eq_constraints, solver.linear_int_constraints)
#     linear_constraints_evaluation = evaluateScalarAffineFunctionConstraints(linear_constraints, x)
#
#     quadratic_constraints = vcat(solver.quadratic_le_constraints, solver.quadratic_ge_constraints, solver.quadratic_eq_constraints, solver.quadratic_int_constraints)
#     quadratic_constraints_evaluation = evaluateScalarQuadraticFunctionConstraints(quadratic_constraints, x)
#
#     non_linear_constraints_evaluation = zeros(Float64, length(m.solver.nlp_data.evaluator.constraints))
#     MOI.eval_constraint(m.solver.nlp_data.evaluator, non_linear_constraints_evaluation, x)
#     println("m.solver.nlp_data.evaluator.constraints: ", m.solver.nlp_data.evaluator.constraints)
#     println("m.solver.nlp_data.evaluator.constraints: ", length(m.solver.nlp_data.evaluator.constraints))
#     a = vcat(linear_constraints_evaluation, quadratic_constraints_evaluation, non_linear_constraints_evaluation)
#     #a = cons(m.nlp, _cute_x(m, x) )
#     if a == Float64[]
#       a = zeros(1)
#     end
#     return [lb(a, m.bcon); ub(a, m.bcon); lb(x, m.bvar); ub(x, m.bvar)];
# end

function _cute_x(m::Class_CUTEst, x::Array{Float64,1})
    if m.nlp != nothing
        ind = _i_not_fixed(m.nlp)
        #@show length(ind)
        if(length(x) != length(ind))
          error("$(length(x)) = length(x) != length(i_not_fixed) = $(length(ind))")
        end


        cute_x = deepcopy(m.nlp.meta.lvar) # get correct values of fixed variables
        cute_x[ind] = x

        return cute_x
    end

    ind = _i_not_fixed(m.solver.variable_info)
    #@show length(ind)
    if(length(x) != length(ind))
      error("$(length(x)) = length(x) != length(i_not_fixed) = $(length(ind))")
    end

    indexes = [i for i in 1:length(m.solver.variable_info)]
    lvar, uvar = extract_lvar_uvar(m.solver.variable_info, indexes)
    cute_x = deepcopy(lvar) # get correct values of fixed variables
    cute_x[ind] = x

    return cute_x
end

#IMPORTANT_
# function _cute_x(m::Class_CUTEst, x::Array{Float64,1})
#     if m.nlp != nothing
#         ind = _i_not_fixed(m.nlp)
#         #@show length(ind)
#         if(length(x) != length(ind))
#           error("$(length(x)) = length(x) != length(i_not_fixed) = $(length(ind))")
#         end
#
#
#         cute_x = deepcopy(m.nlp.meta.lvar) # get correct values of fixed variables
#         cute_x[ind] = x
#
#         return cute_x
#     end
#
#     ind = _i_not_fixed(m.solver.variable_info)
#     #@show length(ind)
#     if(length(x) != length(ind))
#       error("$(length(x)) = length(x) != length(i_not_fixed) = $(length(ind))")
#     end
#
#     indexes = [i for i in 1:length(m.solver.variable_info)]
#     lvar, uvar = extract_lvar_uvar(m.solver.variable_info, indexes)
#     println("lvar: ", lvar)
#     cute_x = deepcopy(lvar) # get correct values of fixed variables
#     cute_x[ind] = x
#
#     return cute_x
# end

function eval_jac(m::Class_CUTEst, x::Array{Float64,1})
    if m.nlp != nothing
        cute_x = _cute_x(m, x)
        J_full_T = jac(m.nlp, cute_x)'
        J_T = J_full_T[_i_not_fixed(m.nlp),:];
        #my_eye = SparseArrays.speye(length(x))
        my_eye = SparseMatrixCSC{Float64}(LinearAlgebra.I, length(x), length(x))
        Q_T = [J_T[:,m.bcon.l_i] -J_T[:,m.bcon.u_i] my_eye[:,m.bvar.l_i] -my_eye[:,m.bvar.u_i]];
        return SparseArrays.sparse(Q_T')
    end

    cute_x = _cute_x(m, x)
    solver = m.solver
    linear_constraints = vcat(solver.linear_le_constraints, solver.linear_ge_constraints, solver.linear_eq_constraints, solver.linear_int_constraints)
    linear_constraints_jac_structure = evaluateScalarAffineFunctionConstraintsJacobianStructure(linear_constraints, x)

    quadratic_constraints = vcat(solver.quadratic_le_constraints, solver.quadratic_ge_constraints, solver.quadratic_eq_constraints, solver.quadratic_int_constraints)
    quadratic_constraints_jac_structure = evaluateScalarQuadraticFunctionConstraintsJacobianStructure(quadratic_constraints, x)
    d = m.solver.nlp_data.evaluator
    jac_struct = MOI.jacobian_structure(d)
    nnz_jac = length(jac_struct)
    J = zeros(Float64, nnz_jac)
    MOI.eval_constraint_jacobian(m.solver.nlp_data.evaluator, J, cute_x)

    J_full_T = nothing
    if linear_constraints_jac_structure!= nothing && !isempty(linear_constraints_jac_structure)
        J_full_T = hcat(linear_constraints_jac_structure)
        if quadratic_constraints_jac_structure != nothing && !isempty(quadratic_constraints_jac_structure)
            J_full_T = hcat(J_full_T, quadratic_constraints_jac_structure)
        end
        if J != nothing && !isempty(J)
            J_full_T = hcat(J_full_T, J)
        end
    elseif quadratic_constraints_jac_structure != nothing && !isempty(quadratic_constraints_jac_structure)
        J_full_T = hcat(quadratic_constraints_jac_structure)
        if J != nothing && !isempty(J)
            J_full_T = hcat(J_full_T, J)
        end
    elseif J != nothing && !isempty(J)
        count_variables = length(m.solver.variable_info)
        count_nonlinear_constraints = number_constraints_non_linear(m.solver)
        J_full_T = zeros(count_nonlinear_constraints, count_variables)
        for i in 1:length(jac_struct)
            J_full_T[jac_struct[i][1], jac_struct[i][2]] = J[i]
        end
        J_full_T = J_full_T'
    end

    J_T = J_full_T[_i_not_fixed(m.solver.variable_info),:];
    my_eye = SparseMatrixCSC{Float64}(LinearAlgebra.I, length(x), length(x))
    Q_T = [J_T[:,m.bcon.l_i] -J_T[:,m.bcon.u_i] my_eye[:,m.bvar.l_i] -my_eye[:,m.bvar.u_i]];
    return SparseArrays.sparse(Q_T')
end

#IMPORTANT_
# function eval_jac(m::Class_CUTEst, x::Array{Float64,1})
#     if m.nlp != nothing
#         cute_x = _cute_x(m, x)
#         J_full_T = jac(m.nlp, cute_x)'
#         J_T = J_full_T[_i_not_fixed(m.nlp),:];
#         #my_eye = SparseArrays.speye(length(x))
#         my_eye = SparseMatrixCSC{Float64}(LinearAlgebra.I, length(x), length(x))
#         Q_T = [J_T[:,m.bcon.l_i] -J_T[:,m.bcon.u_i] my_eye[:,m.bvar.l_i] -my_eye[:,m.bvar.u_i]];
#         return SparseArrays.sparse(Q_T')
#     end
#
#     cute_x = _cute_x(m, x)
#     solver = m.solver
#     linear_constraints = vcat(solver.linear_le_constraints, solver.linear_ge_constraints, solver.linear_eq_constraints, solver.linear_int_constraints)
#     linear_constraints_jac_structure = evaluateScalarAffineFunctionConstraintsJacobianStructure(linear_constraints, x)
#
#     quadratic_constraints = vcat(solver.quadratic_le_constraints, solver.quadratic_ge_constraints, solver.quadratic_eq_constraints, solver.quadratic_int_constraints)
#     quadratic_constraints_jac_structure = evaluateScalarQuadraticFunctionConstraintsJacobianStructure(quadratic_constraints, x)
#     d = m.solver.nlp_data.evaluator
#     jac_struct = MOI.jacobian_structure(d)
#     nnz_jac = length(jac_struct)
#     J = zeros(Float64, nnz_jac)
#     MOI.eval_constraint_jacobian(m.solver.nlp_data.evaluator, J, cute_x)
#
#     println("linear_constraints_jac_structure: ", linear_constraints_jac_structure)
#     println("quadratic_constraints_jac_structure: ", quadratic_constraints_jac_structure)
#     println("J: ", J)
#
#     J_full_T = nothing
#     if !isempty(linear_constraints_jac_structure)
#         J_full_T = hcat(linear_constraints_jac_structure)
#         if !isempty(quadratic_constraints_jac_structure)
#             J_full_T = hcat(J_full_T, quadratic_constraints_jac_structure)
#         end
#         if !isempty(J)
#             J_full_T = hcat(J_full_T, J)
#         end
#     elseif !isempty(quadratic_constraints_jac_structure)
#         J_full_T = hcat(quadratic_constraints_jac_structure)
#         if !isempty(J)
#             J_full_T = hcat(J_full_T, J)
#         end
#     elseif !isempty(J)
#         J_full_T = hcat(J)
#     end
#
#     # J_full_T = J'
#     J_T = J_full_T[_i_not_fixed(m.solver.variable_info),:];
#     println("J_T: ", J_T)
#     #my_eye = SparseArrays.speye(length(x))
#     my_eye = SparseMatrixCSC{Float64}(LinearAlgebra.I, length(x), length(x))
#     Q_T = [J_T[:,m.bcon.l_i] -J_T[:,m.bcon.u_i] my_eye[:,m.bvar.l_i] -my_eye[:,m.bvar.u_i]];
#     return SparseArrays.sparse(Q_T')
#
#     #return @time [J[m.bcon.l_i,:]; -J[m.bcon.u_i,:]; my_eye[m.bvar.l_i,:]; -my_eye[m.bvar.u_i,:]];
# end

function eval_grad_f(m::Class_CUTEst, x::Array{Float64,1})
    if m.nlp != nothing
        return grad(m.nlp, _cute_x(m, x))[_i_not_fixed(m.nlp)]
    end

    g = zeros(length(x))
    MOI.eval_objective_gradient(m.solver.nlp_data.evaluator, g, _cute_x(m, x))
    return g[_i_not_fixed(m.solver.variable_info)]
end

function get_reducedcosts(m::Class_CUTEst, y::Array{Float64,1})
    if m.nlp != nothing
        rc = zeros(m.nlp.meta.nvar)
        st_l = ncons_orginal(m)
        st_u = length(m.bvar.l) + ncons_orginal(m)
        rc[m.bvar.l_i] += y[(st_l+1):st_u]
        rc[m.bvar.u_i] -= y[(st_u+1):(length(m.bvar.u) + st_u)]
        return rc
    end

    rc = zeros(length(m.solver.variable_info))
    st_l = ncons_orginal(m)
    st_u = length(m.bvar.l) + ncons_orginal(m)
    rc[m.bvar.l_i] += y[(st_l+1):st_u]
    rc[m.bvar.u_i] -= y[(st_u+1):(length(m.bvar.u) + st_u)]
    return rc
end

function y_cons_net(m::Class_CUTEst, y::Array{Float64,1})
    return;
end

# function eval_jtprod(m::Class_CUTEst, x::Array{Float64,1}, y::Array{Float64,1})
#     return jtprod(nlp_raw, x, y)
# end

function eval_schur_and_jtprod(m::Class_CUTEst, x::Array{Float64,1}, y::Array{Float64,1}, s::Array{Float64,1})

end

#=function ncon(m::Class_CUTEst)
    return length(m.bcon.l) + length(m.bcon.u)
end=#

#=function nvar(m::Class_CUTEst)
    return m.nlp.meta.nvar
end=#
#=
function make_symmetric(M::SparseMatrixCSC{Float64,Int32})
    n = size(M,1)
    new_M = spzeros(n,n)
    rows = rowvals(M)
    vals = nonzeros(M)
    for col = 1:n
      for j in nzrange(M, col)
         row = rows[j]
         val = vals[j]
         # perform sparse wizardry...
         new_M[col,row] = vals[j]
         new_M[row,col] = vals[j]
      end
    end

    return new_M
end=#

function eval_lag_hess(m::Class_CUTEst, x::Array{Float64,1}, y::Array{Float64,1}, w::Float64)
    if m.nlp != nothing
        y_cons = zeros(m.nlp.meta.ncon)
        y_cons[m.bcon.l_i] -= y_l_con(y, m)
        y_cons[m.bcon.u_i] += y_u_con(y, m)
        H = hess(m.nlp, _cute_x(m, x), y_cons; obj_weight=w);

        H = LowerTriangular(H)
        ind = _i_not_fixed(m.nlp)
        H_not_fixed = H[ind,ind]
        return H_not_fixed
    end

    y_cons = zeros(number_constraints(m.solver))
    y_cons[m.bcon.l_i] -= y_l_con(y, m)
    y_cons[m.bcon.u_i] += y_u_con(y, m)

    d = m.solver.nlp_data.evaluator
    hessian_sparsity = MOI.hessian_lagrangian_structure(d)
    #Sparse hessian entry vector
    H_vec = zeros(Float64, length(hessian_sparsity))
    MOI.eval_hessian_lagrangian(d, H_vec, _cute_x(m, x), w, y_cons)
    H = spzeros(length(m.solver.variable_info), length(m.solver.variable_info))
    if !isempty(hessian_sparsity)
        index = 1
        if length(hessian_sparsity) < 3
            index = length(hessian_sparsity)
        elseif length(m.solver.variable_info) > 1
            index = 3
            for i in 3:length(m.solver.variable_info)
                index = index + i
            end
        end
        for i in 1:index
            H[hessian_sparsity[i][1], hessian_sparsity[i][2]]= H_vec[i]
        end
    end
    ind = _i_not_fixed(m.solver.variable_info)
    H_not_fixed = H[ind,ind]
    return H_not_fixed
end

#IMPORTANT_
# function eval_lag_hess(m::Class_CUTEst, x::Array{Float64,1}, y::Array{Float64,1}, w::Float64)
#     y_cons = zeros(number_constraints(m.solver))
#     println("_____y: ", y)
#     println("_____y_l_con(y, m)", y_l_con(y, m))
#     y_cons[m.bcon.l_i] -= y_l_con(y, m)
#     y_cons[m.bcon.u_i] += y_u_con(y, m)
#
#     #H = hess(m.nlp, _cute_x(m, x), obj_weight=w, y=y_cons);
#     #H = hess(m.nlp, _cute_x(m, x), y=y_cons; obj_weight=w);
#     #H = hess(m.nlp, _cute_x(m, x), y_cons; obj_weight=w);
#     #After NLPModels 0.16.0, they started returnning Symmetric Hessian instead of the Lower Triangular
#     #H = LowerTriangular(H)
#     # println("_____________________________", H)
#     # println("_____________________________", LowerTriangular(H))
#     #H = hess(m.nlp, _cute_x(m, x), w, y_cons);
#     #@show typeof(H)
#     H = spzeros(length(m.solver.variable_info), length(m.solver.variable_info))
#     MOI.eval_hessian_lagrangian(m.solver.nlp_data.evaluator, H, _cute_x(m, x), w, y_cons)
#
#     ind = _i_not_fixed(m.solver.variable_info)
#
#     H_not_fixed = H[ind,ind]
#     #H_true = SparseArrays.sparse(Symmetric(H_not_fixed,:L) + spzeros(length(x),length(x)))
#     #H_true = H_not_fixed + H_not_fixed' - spdiagm(diag(H_not_fixed))
#     #convert()
#     #@time H_true = make_symmetric(H_not_fixed)
#     #@show LinearAlgebra.norm(H_true - H2,Inf)
#
#     return H_not_fixed
# end

# function eval_Jt_prod(m::Class_CUTEst, x::Array{Float64,1}, y::Array{Float64,1})
#     return jtprod(m.nlp, x, y)
# end
