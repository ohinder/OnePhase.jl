#using MathProgBase
import MathOptInterface
#using SparseArrays

########################################################
## this code is based on ModelReader in NLPModels
## and KNITRO.jl
########################################################

export OnePhaseSolver

const MOI  = MathOptInterface
const MOIU = MathOptInterface.Utilities

# TODO
const SF = Union{MOI.ScalarAffineFunction{Float64},
                 MOI.ScalarQuadraticFunction{Float64}}
const VAF = MOI.VectorAffineFunction{Float64}
const VOV = MOI.VectorOfVariables

# ScalarAffineFunctions and VectorAffineFunctions
const SAF = MOI.ScalarAffineFunction{Float64}
const AF = Union{SAF, VAF}

const SS = Union{MOI.EqualTo{Float64},
                 MOI.GreaterThan{Float64},
                 MOI.LessThan{Float64},
                 MOI.Interval{Float64}}
# LinSets
const LS = Union{MOI.EqualTo{Float64},
                 MOI.GreaterThan{Float64},
                 MOI.LessThan{Float64}}
# VecLinSets
const VLS = Union{MOI.Nonnegatives,
                  MOI.Nonpositives,
                  MOI.Zeros}
##################################################
mutable struct VariableInfo
    lower_bound::Float64  # May be -Inf even if has_lower_bound == true
    has_lower_bound::Bool # Implies lower_bound == Inf
    lower_bound_dual_start::Union{Nothing, Float64}
    upper_bound::Float64  # May be Inf even if has_upper_bound == true
    has_upper_bound::Bool # Implies upper_bound == Inf
    upper_bound_dual_start::Union{Nothing, Float64}
    is_fixed::Bool        # Implies lower_bound == upper_bound and !has_lower_bound and !has_upper_bound.
    start::Union{Nothing, Float64}
end
VariableInfo() = VariableInfo(-Inf, false, nothing, Inf, false, nothing, false, nothing)

mutable struct OnePhaseProblem
    status::Symbol  # Final status
    #status::MOI.TerminationStatusCode
    x::Vector{Float64}  # Starting and final solution
    lambda::Vector{Float64}
    g::Vector{Float64}  # Final constraint values
    obj_val::Float64  # (length 1) Final objective
    solve_time::Float64

    # Custom attributes of the OnePhaseSolver
    iter::Class_iterate
    hist::Array{alg_history2,1}
    pars::Class_parameters

    function OnePhaseProblem()
        return new()
    end
end

##################################################
# EmptyNLPEvaluator for non-NLP problems.
struct EmptyNLPEvaluator <: MOI.AbstractNLPEvaluator end
MOI.features_available(::EmptyNLPEvaluator) = [:Grad, :Jac, :Hess]
MOI.initialize(::EmptyNLPEvaluator, features) = nothing
MOI.eval_objective(::EmptyNLPEvaluator, x) = NaN
function MOI.eval_constraint(::EmptyNLPEvaluator, g, x)
    @assert length(g) == 0
    return
end

MOI.eval_objective_gradient(::EmptyNLPEvaluator, g, x) = nothing
MOI.jacobian_structure(::EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
MOI.hessian_lagrangian_structure(::EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
function MOI.eval_constraint_jacobian(::EmptyNLPEvaluator, J, x)
    @assert length(J) == 0
    return
end
function MOI.eval_hessian_lagrangian(::EmptyNLPEvaluator, H, x, s, mu)
    @assert length(H) == 0
    return
end
function MOI.eval_hessian_lagrangian(::Nothing, H, x, s, mu)
    @assert length(H) == 0
    return
end

empty_nlp_data() = MOI.NLPBlockData([], EmptyNLPEvaluator(), false)

mutable struct ConstraintInfo{F, S}
    func::F
    set::S
    dual_start::Union{Nothing, Float64}
end

ConstraintInfo(func, set) = ConstraintInfo(func, set, nothing)

mutable struct OnePhaseSolver <: MOI.AbstractOptimizer
    #inner::OnePhaseProblem
	inner::Union{OnePhaseProblem, Nothing}

    # Problem data.
	variable_info::Vector{VariableInfo}
	nlp_data::MOI.NLPBlockData
	sense :: MOI.OptimizationSense
	objective::Union{MOI.VariableIndex, MOI.ScalarAffineFunction{Float64}, MOI.ScalarQuadraticFunction{Float64}, Nothing}
    linear_le_constraints::Vector{ConstraintInfo{MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}}}
    linear_ge_constraints::Vector{ConstraintInfo{MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}}}
    linear_eq_constraints::Vector{ConstraintInfo{MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}}}
	linear_int_constraints::Vector{ConstraintInfo{MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64}}}
    quadratic_le_constraints::Vector{ConstraintInfo{MOI.ScalarQuadraticFunction{Float64}, MOI.LessThan{Float64}}}
    quadratic_ge_constraints::Vector{ConstraintInfo{MOI.ScalarQuadraticFunction{Float64}, MOI.GreaterThan{Float64}}}
    quadratic_eq_constraints::Vector{ConstraintInfo{MOI.ScalarQuadraticFunction{Float64}, MOI.EqualTo{Float64}}}
	quadratic_int_constraints::Vector{ConstraintInfo{MOI.ScalarQuadraticFunction{Float64}, MOI.Interval{Float64}}}
    nlp_dual_start::Union{Nothing, Vector{Float64}}

	##Constraint mappings
	constraint_mapping::Dict{MOI.ConstraintIndex, Union{Cint, Vector{Cint}}}
	# Parameters.
	silent::Bool
	options::Dict{String, Any}

	# Solution attributes.
    solve_time::Float64
end

function OnePhaseSolver(; options...)
	options_dict = Dict{String, Any}()

	for (name, value) in options
        options_dict[string(name)] = value
    end

	onePhaseSolverModel = OnePhaseSolver(
        OnePhaseProblem(),
        [],
        empty_nlp_data(),
        MOI.FEASIBILITY_SENSE,
        nothing,
        [],
        [],
        [],
        [],
        [],
        [],
		[],
		[],
        nothing,
		Dict{MOI.ConstraintIndex, Int}(),
        false,
        options_dict,
        NaN,
    )
    set_options(onePhaseSolverModel, options)

	return onePhaseSolverModel
end

function set_options(model::OnePhaseSolver, options)
    for (name, value) in options
        sname = string(name)
        MOI.set(model, MOI.RawOptimizerAttribute(sname), value)
    end
    return
end

MOI.get(::OnePhaseSolver, ::MOI.SolverName) = "OnePhaseSolver"

"""
    MOI.is_empty(model::OnePhaseSolver )
"""

function MOI.is_empty(model::OnePhaseSolver)
    return isempty(model.variable_info) &&
           model.nlp_data.evaluator isa EmptyNLPEvaluator &&
           model.sense == MOI.FEASIBILITY_SENSE &&
           isempty(model.linear_le_constraints) &&
           isempty(model.linear_ge_constraints) &&
           isempty(model.linear_eq_constraints) &&
		   isempty(model.linear_int_constraints) &&
           isempty(model.quadratic_le_constraints) &&
           isempty(model.quadratic_ge_constraints) &&
           isempty(model.quadratic_eq_constraints) &&
		   isempty(model.quadratic_int_constraints)
end

# MOI.get(model::OnePhaseSolver, ::MOI.SolveTime) = model.solve_time

function MOI.empty!(model::OnePhaseSolver)
    model.inner = nothing
    empty!(model.variable_info)
    model.nlp_data = empty_nlp_data()
    model.sense = MOI.FEASIBILITY_SENSE
    model.objective = nothing
    empty!(model.linear_le_constraints)
    empty!(model.linear_ge_constraints)
    empty!(model.linear_eq_constraints)
	empty!(model.linear_int_constraints)
    empty!(model.quadratic_le_constraints)
    empty!(model.quadratic_ge_constraints)
    empty!(model.quadratic_eq_constraints)
	empty!(model.quadratic_int_constraints)
    model.nlp_dual_start = nothing
end

function MOI.get(model::OnePhaseSolver, ::MOI.ListOfVariableIndices)
	return [MOI.VariableIndex(i) for i in 1:length(model.variable_info)]
end

function MOI.add_variable(model::OnePhaseSolver)
    push!(model.variable_info, VariableInfo())
    return MOI.VariableIndex(length(model.variable_info))
end

function MOI.add_variables(model::OnePhaseSolver, n::Int)
    return [MOI.add_variable(model) for i in 1:n]
end

function MOI.is_valid(model::OnePhaseSolver, vi::MOI.VariableIndex)
    return vi.value in eachindex(model.variable_info)
end

function MOI.Utilities.supports_default_copy_to(::OnePhaseSolver, copy_names::Bool)
    return !copy_names
end

function MOI.copy_to(model::OnePhaseSolver, src::MOI.ModelLike; copy_names = false)
    return MOI.Utilities.default_copy_to(model, src, copy_names)
end
# MathOptInterface constraints

##################################################
## Support constraints
MOI.supports_constraint(::OnePhaseSolver, ::Type{MOI.VariableIndex}, ::Type{<:SS}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{MOI.LessThan{Float64}}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{MOI.GreaterThan{Float64}}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{MOI.EqualTo{Float64}}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{MOI.Interval{Float64}}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{MOI.ScalarQuadraticFunction{Float64}}, ::Type{MOI.LessThan{Float64}}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{MOI.ScalarQuadraticFunction{Float64}}, ::Type{MOI.GreaterThan{Float64}}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{MOI.ScalarQuadraticFunction{Float64}}, ::Type{MOI.EqualTo{Float64}}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{MOI.ScalarQuadraticFunction{Float64}}, ::Type{MOI.Interval{Float64}}) = true

MOI.supports_constraint(::OnePhaseSolver, ::Type{<:SF}, ::Type{<:SS}) = true

function has_upper_bound(model::OnePhaseSolver, vi::MOI.VariableIndex)
    return model.variable_info[vi.value].has_upper_bound
end

function has_lower_bound(model::OnePhaseSolver, vi::MOI.VariableIndex)
    return model.variable_info[vi.value].has_lower_bound
end

function is_fixed(model::OnePhaseSolver, vi::MOI.VariableIndex)
    return model.variable_info[vi.value].is_fixed
end

function MOI.add_constraint(
    model::OnePhaseSolver, v::MOI.VariableIndex, lt::MOI.LessThan{Float64},
)
	vi = v
    MOI.throw_if_not_valid(model, vi)
    if isnan(lt.upper)
        error("Invalid upper bound value $(lt.upper).")
    end
    if has_upper_bound(model, vi)
        throw(MOI.UpperBoundAlreadySet{typeof(lt), typeof(lt)}(vi))
    end
    if is_fixed(model, vi)
        throw(MOI.UpperBoundAlreadySet{MOI.EqualTo{Float64}, typeof(lt)}(vi))
    end

	col = vi.value
    model.variable_info[col].upper_bound = lt.upper
    model.variable_info[col].has_upper_bound = true
	ci = MOI.ConstraintIndex{MOI.VariableIndex, MOI.LessThan{Float64}}(col)
    return ci
end

function MOI.set(
    model::OnePhaseSolver,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VariableIndex, MOI.LessThan{Float64}},
    set::MOI.LessThan{Float64},
)
    MOI.throw_if_not_valid(model, ci)
    model.variable_info[ci.value].upper_bound = set.upper
    return
end

function MOI.delete(
    model::OnePhaseSolver,
    ci::MOI.ConstraintIndex{MOI.VariableIndex, MOI.LessThan{Float64}},
)
    MOI.throw_if_not_valid(model, ci)
    model.variable_info[ci.value].upper_bound = Inf
    model.variable_info[ci.value].has_upper_bound = false
    return
end

function MOI.add_constraint(
    model::OnePhaseSolver, v::MOI.VariableIndex, gt::MOI.GreaterThan{Float64},
)
	vi = v
    MOI.throw_if_not_valid(model, vi)
    if isnan(gt.lower)
        error("Invalid lower bound value $(gt.lower).")
    end
    if has_lower_bound(model, vi)
        throw(MOI.LowerBoundAlreadySet{typeof(gt), typeof(gt)}(vi))
    end
    if is_fixed(model, vi)
        throw(MOI.LowerBoundAlreadySet{MOI.EqualTo{Float64}, typeof(gt)}(vi))
    end

	col = vi.value
    model.variable_info[col].lower_bound = gt.lower
    model.variable_info[col].has_lower_bound = true
	ci = MOI.ConstraintIndex{MOI.VariableIndex, MOI.GreaterThan{Float64}}(col)
	return ci
end

function MOI.set(
    model::OnePhaseSolver,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VariableIndex, MOI.GreaterThan{Float64}},
    set::MOI.GreaterThan{Float64},
)
    MOI.throw_if_not_valid(model, ci)
    model.variable_info[ci.value].lower_bound = set.lower
    return
end

function MOI.delete(
    model::OnePhaseSolver,
    ci::MOI.ConstraintIndex{MOI.VariableIndex, MOI.GreaterThan{Float64}},
)
    MOI.throw_if_not_valid(model, ci)
    model.variable_info[ci.value].lower_bound = -Inf
    model.variable_info[ci.value].has_lower_bound = false
    return
end

function MOI.add_constraint(
    model::OnePhaseSolver,  v::MOI.VariableIndex, eq::MOI.EqualTo{Float64},
)
	vi = v
    MOI.throw_if_not_valid(model, vi)
    if isnan(eq.value)
        error("Invalid fixed value $(eq.value).")
    end
    if has_lower_bound(model, vi)
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{Float64}, typeof(eq)}(vi))
    end
    if has_upper_bound(model, vi)
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{Float64}, typeof(eq)}(vi))
    end
    if is_fixed(model, vi)
        throw(MOI.LowerBoundAlreadySet{typeof(eq), typeof(eq)}(vi))
    end

	col = vi.value
    model.variable_info[col].lower_bound = eq.value
    model.variable_info[col].upper_bound = eq.value
    model.variable_info[col].is_fixed = true
	ci = MOI.ConstraintIndex{MOI.VariableIndex, MOI.EqualTo{Float64}}(col)
	return ci
end

function MOI.set(
    model::OnePhaseSolver,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VariableIndex, MOI.EqualTo{Float64}},
    set::MOI.EqualTo{Float64},
)
    MOI.throw_if_not_valid(model, ci)
    model.variable_info[ci.value].lower_bound = set.value
    model.variable_info[ci.value].upper_bound = set.value
    return
end

function MOI.delete(
    model::OnePhaseSolver,
    ci::MOI.ConstraintIndex{MOI.VariableIndex, MOI.EqualTo{Float64}},
)
    MOI.throw_if_not_valid(model, ci)
    model.variable_info[ci.value].lower_bound = -Inf
    model.variable_info[ci.value].upper_bound = Inf
    model.variable_info[ci.value].is_fixed = false
    return
end

macro define_add_constraint(function_type, set_type, prefix)
    array_name = Symbol(string(prefix) * "_constraints")
    return quote
        function MOI.add_constraint(
            model::OnePhaseSolver, func::$function_type, set::$set_type,
        )
            check_inbounds(model, func)
            push!(model.$(array_name), ConstraintInfo(func, set))

			col = length(model.$(array_name))
			ci = MOI.ConstraintIndex{$function_type, $set_type}(col)
			model.constraint_mapping[ci] = convert(Cint, col)
			return ci
        end
    end
end

@define_add_constraint(
    MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}, linear_le,
)

@define_add_constraint(
    MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}, linear_ge,
)

@define_add_constraint(
    MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}, linear_eq,
)

@define_add_constraint(
    MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64}, linear_int,
)

@define_add_constraint(
    MOI.ScalarQuadraticFunction{Float64}, MOI.LessThan{Float64}, quadratic_le,
)

@define_add_constraint(
    MOI.ScalarQuadraticFunction{Float64}, MOI.GreaterThan{Float64}, quadratic_ge,
)

@define_add_constraint(
    MOI.ScalarQuadraticFunction{Float64}, MOI.EqualTo{Float64}, quadratic_eq,
)

@define_add_constraint(
    MOI.ScalarQuadraticFunction{Float64}, MOI.Interval{Float64}, quadratic_int,
)

function MOI.set(
    model::OnePhaseSolver,
    ::MOI.ObjectiveFunction,
    func::Union{
        MOI.VariableIndex,
        MOI.ScalarAffineFunction,
        MOI.ScalarQuadraticFunction,
    },
)
    check_inbounds(model, func)
    model.objective = func
    return
end

function NLPModels.obj(nlp::MathOptNLPModel, x::AbstractVector)
  NLPModels.increment!(nlp, :neval_obj)
  if nlp.obj.type == "LINEAR"
    res = dot(nlp.obj.gradient, x) + nlp.obj.constant
  end
  if nlp.obj.type == "QUADRATIC"
    res =
      0.5 * coo_sym_dot(nlp.obj.hessian.rows, nlp.obj.hessian.cols, nlp.obj.hessian.vals, x, x) +
      dot(nlp.obj.gradient, x) +
      nlp.obj.constant
  end
  if nlp.obj.type == "NONLINEAR"
    res = MOI.eval_objective(nlp.eval, x)
  end
  return res
end

function append_to_hessian_sparsity!(
    ::Any,
    ::Union{MOI.VariableIndex,MOI.ScalarAffineFunction},
)
    return nothing
end

function append_to_hessian_sparsity!(
    hessian_sparsity,
    quad::MOI.ScalarQuadraticFunction,
)
    for term in quad.quadratic_terms
        push!(
            hessian_sparsity,
            # (term.variable_index_1.value, term.variable_index_2.value),
			(term.variable_1.value, term.variable_2.value),
        )
    end
end

function hessian_lagrangian_structure(model::OnePhaseSolver, nlp :: MathOptNLPModel)
    hessian_sparsity = Tuple{Int64,Int64}[]
    if model.objective !== nothing
        append_to_hessian_sparsity!(hessian_sparsity, model.obj)
    end
    for info in model.quadratic_le_constraints
        append_to_hessian_sparsity!(hessian_sparsity, info.func)
    end
    for info in model.quadratic_ge_constraints
        append_to_hessian_sparsity!(hessian_sparsity, info.func)
    end
    for info in model.quadratic_eq_constraints
        append_to_hessian_sparsity!(hessian_sparsity, info.func)
    end
    nlp_hessian_sparsity =
        MOI.hessian_lagrangian_structure(nlp.eval)
    append!(hessian_sparsity, nlp_hessian_sparsity)
    return hessian_sparsity
end

function hess_coord(nlp :: MathOptNLPModel, x :: Array{Float64};
  obj_weight :: Float64=1.0, y :: Array{Float64}=zeros(nlp.meta.ncon))
  NLPModels.increment!(nlp, :neval_hess)
  MOI.eval_hessian_lagrangian(nlp.eval, nlp.obj.hessian.vals, x, obj_weight, y)

  return (NLPModels.hess_structure(nlp)[1], NLPModels.hess_structure(nlp)[2], NLPModels.hess_coord(nlp, x, y, obj_weight=obj_weight))
end

function NLPModels.hess_structure!(
  nlp::MathOptNLPModel,
  rows::AbstractVector{<:Integer},
  cols::AbstractVector{<:Integer},
)
  if nlp.obj.type == "QUADRATIC"
    for index = 1:(nlp.obj.nnzh)
      rows[index] = nlp.obj.hessian.rows[index]
      cols[index] = nlp.obj.hessian.cols[index]
    end
  end
  if (nlp.obj.type == "NONLINEAR") || (nlp.meta.nnln > 0)
    hesslag_struct = MOI.hessian_lagrangian_structure(nlp.eval)
    for index = (nlp.obj.nnzh + 1):(nlp.meta.nnzh)
      shift_index = index - nlp.obj.nnzh
      rows[index] = hesslag_struct[shift_index][1]
      cols[index] = hesslag_struct[shift_index][2]
    end
  end
  return rows, cols
end

############################
## END ModelReader CODE
############################

function status_One_Phase_To_JuMP(status::Symbol)
    # since our status are not equal to JuMPs we need to do a conversion
    if status == :Optimal
        return :Optimal
    elseif status == :primal_infeasible
        return :Infeasible
    elseif status == :dual_infeasible
        return :Unbounded
    elseif status == :MAX_IT || status === :MAX_TIME
        return :UserLimit
    else
        return :Error
    end
end

function create_pars_JuMP(options )
    pars = Class_parameters()
    for (param,value) in options
      what = split(String(param),"!") # we represent a parameter such as init.mu_scale as init!mu_scale because we cannot pass init.mu_scale as a parameter
      node = pars # root
      for i in 1:length(what)
          field = what[i]
          if i < length(what)
              node = getfield(node,Symbol(field))
          else # at the leaf
              setfield!(node,Symbol(field),value)
          end
      end
    end

    return pars
end

MOI.get(model::OnePhaseSolver, ::MOI.RawSolver) = model

# copy
MOI.supports_incremental_interface(solver::OnePhaseSolver) = true
function MOI.copy_to(solver::OnePhaseSolver, src::MOI.ModelLike)
    return MOI.Utilities.automatic_copy_to(model, src)
end

function MOI.optimize!(solver :: OnePhaseSolver)
    t = time()

	has_nlp_objective = false
	if !isa(solver.nlp_data.evaluator, EmptyNLPEvaluator)
		features = MOI.features_available(solver.nlp_data.evaluator)::Vector{Symbol}
		has_hessian = (:Hess in features)
		has_hessvec = (:HessVec in features)
		has_nlp_objective = solver.nlp_data.has_objective
		num_nlp_constraints = length(solver.nlp_data.constraint_bounds)
		has_nlp_constraints = (num_nlp_constraints > 0)

		init_feat = Symbol[]
		has_nlp_objective && push!(init_feat, :Grad)
		if has_hessian
			push!(init_feat, :Hess)
		elseif has_hessvec
			push!(init_feat, :HessVec)
		end

		if has_nlp_constraints
			push!(init_feat, :Jac)
		end

		MOI.initialize(solver.nlp_data.evaluator, init_feat)
	end

    pars = create_pars_JuMP(solver.options)

    iter, status, hist, t, err, timer = one_phase_solve(solver, pars)
    solver.inner = OnePhaseProblem()
    solver.inner.status = status_One_Phase_To_JuMP(status)
    solver.inner.x = get_original_x(iter)
    solver.inner.obj_val = iter.cache.fval
    solver.inner.lambda = get_y(iter)
    solver.inner.solve_time = time() - t

    # custom one phase features
    solver.inner.pars = pars
    solver.inner.iter = iter
    solver.inner.hist = hist
end

function MOI.optimize!(solver :: OnePhaseSolver, jumpModel:: Model)
    t = time()
    nlp = MathOptNLPModel(jumpModel)
    pars = create_pars_JuMP(solver.options)

    iter, status, hist, t, err, timer = one_phase_solve(nlp,pars)
    solver.inner = OnePhaseProblem()
    solver.inner.status = status_One_Phase_To_JuMP(status)
    solver.inner.x = get_original_x(iter)
    solver.inner.obj_val = iter.cache.fval
    solver.inner.lambda = get_y(iter)
    solver.inner.solve_time = time() - t

    # custom one phase features
    solver.inner.pars = pars
    solver.inner.iter = iter
    solver.inner.hist = hist
end

function MOI.get(
    model::MathOptInterface.Utilities.CachingOptimizer,
    attr::MOI.TerminationStatus,
)
    return MOI.get(model.optimizer, attr)
end

function check_inbounds(model::OnePhaseSolver, vi::MOI.VariableIndex)
    return MOI.throw_if_not_valid(model, vi)
end

function check_inbounds(model::OnePhaseSolver, aff::MOI.ScalarAffineFunction)
    for term in aff.terms
		MOI.throw_if_not_valid(model, term.variable)
    end
end

function check_inbounds(model::OnePhaseSolver, quad::MOI.ScalarQuadraticFunction)
    for term in quad.affine_terms
		MOI.throw_if_not_valid(model, term.variable)
    end
    for term in quad.quadratic_terms
		MOI.throw_if_not_valid(model, term.variable_1)
        MOI.throw_if_not_valid(model, term.variable_2)
    end
end

MOI.supports(::OnePhaseSolver, ::MOI.NLPBlock) = true

function MOI.supports(
    ::OnePhaseSolver, SF
)
    return true
end

function MOI.supports(
    ::OnePhaseSolver, ::MOI.ObjectiveFunction{MOI.VariableIndex}
)
    return true
end

function MOI.supports(
    ::OnePhaseSolver, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}
)
    return true
end

function MOI.supports(
    ::OnePhaseSolver, ::MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}
)
    return true
end

MOI.supports(::OnePhaseSolver, ::MOI.ObjectiveSense) = true

MOI.supports(::OnePhaseSolver, ::MOI.Silent) = true

MOI.supports(::OnePhaseSolver, ::MOI.RawOptimizerAttribute) = true

function MOI.get(model::OnePhaseSolver, ::MOI.ObjectiveFunction)
    return model.objective
end

function MOI.set(model::OnePhaseSolver, ::MOI.NLPBlock, nlp_data::MOI.NLPBlockData)
    model.nlp_data = nlp_data
    return
end

function MOI.set(
    model::OnePhaseSolver, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense
)
    model.sense = sense
    return
end

MOI.get(model::OnePhaseSolver, ::MOI.ObjectiveSense) = model.sense

function MOI.set(model::OnePhaseSolver, ::MOI.Silent, value)
    model.silent = value
    return
end

MOI.get(model::OnePhaseSolver, ::MOI.Silent) = model.silent

function MOI.supports(
    ::OnePhaseSolver, ::MOI.VariablePrimalStart, ::Type{MOI.VariableIndex}
)
    return true
end

function MOI.set(
    model::OnePhaseSolver,
    ::MOI.VariablePrimalStart,
    vi::MOI.VariableIndex,
    value::Union{Real, Nothing},
)
    MOI.throw_if_not_valid(model, vi)
    model.variable_info[vi.value].start = value
    return
end

function MOI.get(model::MOIU.CachingOptimizer, args...)
	return MOI.get(model.model.optimizer, args)
end

function MOI.get(model::OnePhaseSolver, attr::MOI.VariablePrimal, v::VariableRef)
	return MOI.get(model, attr, v.variable)
end

function MOI.get(
    model::OnePhaseSolver, attr::MOI.VariablePrimal, vi::MOI.VariableIndex,
)
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, vi)
    return model.inner.x[vi.value]
end

function MOI.set(model::OnePhaseSolver, ::MOI.TimeLimitSec, value::Real)
    MOI.set(model, MOI.RawOptimizerAttribute(TIME_LIMIT), Float64(value))
end

function MOI.set(model::OnePhaseSolver, p::MOI.RawOptimizerAttribute, value)
    model.options[p.name] = value
    return
end

function MOI.get(model::OnePhaseSolver, p::MOI.RawOptimizerAttribute)
    if haskey(model.options, p.name)
        return model.options[p.name]
    end
    error("RawParameter with name $(p.name) is not set.")
end

function MOI.get(model::OnePhaseSolver, ::MOI.TerminationStatus)
    if model.inner === nothing
        return MOI.OPTIMIZE_NOT_CALLED
    end
	status = model.inner.status
	return status
end

function MOI.get(model::OnePhaseSolver, ::MOI.RawStatusString)
	return string(model.inner.status)
end

# Ipopt always has an iterate available.
function MOI.get(model::OnePhaseSolver, ::MOI.ResultCount)
    return (model.inner !== nothing) ? 1 : 0
end

function MOI.get(model::OnePhaseSolver, attr::MOI.PrimalStatus)
    if !(1 <= attr.N <= MOI.get(model, MOI.ResultCount()))
        return MOI.NO_SOLUTION
    end

	status = model.inner.status
	return status
end

function MOI.get(model::OnePhaseSolver, attr::MOI.DualStatus)
    if !(1 <= attr.N <= MOI.get(model, MOI.ResultCount()))
        return MOI.NO_SOLUTION
    end
	status = model.inner.status
	return status
end

function MOI.get(model::OnePhaseSolver, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    return model.inner.obj_val
end
