using MathProgBase 
import MathOptInterface
#using SparseArrays

########################################################
## this code is based on ModelReader in NLPModels
## and KNITRO.jl
########################################################

export OnePhaseSolver #, OnePhaseMathProgModel, NonlinearModel

const MOI  = MathOptInterface
const MOIU = MathOptInterface.Utilities

const MPB = MathProgBase

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
    has_lower_bound::Bool # Implies lower_bound == Inf
    has_upper_bound::Bool # Implies upper_bound == Inf
    is_fixed::Bool        # Implies lower_bound == upper_bound and !has_lower_bound and !has_upper_bound.
    name::String
end
VariableInfo() = VariableInfo(false, false, false, "")

mutable struct OnePhaseProblem
    status::Symbol  # Final status

    # For MathProgBase
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
    return
end
MOI.eval_objective_gradient(::EmptyNLPEvaluator, g, x) = nothing
MOI.jacobian_structure(::EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
MOI.hessian_lagrangian_structure(::EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
function MOI.eval_constraint_jacobian(::EmptyNLPEvaluator, J, x)
    return
end
function MOI.eval_hessian_lagrangian(::EmptyNLPEvaluator, H, x, s, mu)
    return
end

empty_nlp_data() = MOI.NLPBlockData([], EmptyNLPEvaluator(), false)


mutable struct OnePhaseSolver <: MOI.AbstractOptimizer
    options::Dict{String, Any}
    #inner::OnePhaseProblem

    #eval :: Union{MOI.AbstractNLPEvaluator, Nothing}
    #numVar::Int
    #numConstr::Int
    #x::Vector{Float64}
    #y::Vector{Float64}
    #lvar::Vector{Float64}
    #uvar::Vector{Float64}
    #lcon::Vector{Float64}
    #ucon::Vector{Float64}
    #sense::Symbol
    #status::Symbol
    nlp_loaded::Bool
    number_solved::Int
end

function OnePhaseSolver(; options...)
    println("#######################", options)
    # Convert Symbol to String in options dictionnary.
    options_dict = Dict{String, Any}()
    for (name, value) in options
        options_dict[string(name)] = value
    end
    onePhaseSolverModel = OnePhaseSolver(Dict(), 0, false)

    set_options(onePhaseSolverModel, options)
    return onePhaseSolverModel
end

function set_options(model::OnePhaseSolver, options)
    for (name, value) in options
        sname = string(name)
        MOI.set(model, MOI.RawParameter(sname), value)
    end
    return
end

#OnePhaseSolver(;kwargs...) = OnePhaseSolver(kwargs)

MOI.get(::OnePhaseSolver, ::MOI.SolverName) = "OnePhaseSolver"

"""
    MOI.is_empty(model::OnePhaseSolver )
"""
#=
function MOI.is_empty(model::OnePhaseSolver)
    #return isempty(model.variable_info) && isempty(model.var_constraints)
	return isempty(model.variable_info) &&
           model.nlp_data.evaluator isa EmptyNLPEvaluator &&
           model.sense == MOI.FEASIBILITY_SENSE &&
           model.number_solved == 0 &&
           isa(model.objective, Nothing) &&
           !model.nlp_loaded
end
=#
function MOI.is_empty(model::OnePhaseSolver)
    return model.number_solved == 0 &&
           #isa(model.objective, Nothing) &&
           !model.nlp_loaded
end


# MathOptInterface constraints

##################################################
## Support constraints

MOI.supports_constraint(::OnePhaseSolver, ::Type{MOI.SingleVariable}, ::Type{MOI.ZeroOne}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{MOI.SingleVariable}, ::Type{MOI.Integer}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{MOI.SingleVariable}, ::Type{<:SS}) = true

MOI.supports_constraint(::OnePhaseSolver, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{MOI.LessThan{Float64}}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{MOI.GreaterThan{Float64}}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{MOI.EqualTo{Float64}}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{MOI.Interval{Float64}}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{MOI.ScalarQuadraticFunction{Float64}}, ::Type{MOI.LessThan{Float64}}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{MOI.ScalarQuadraticFunction{Float64}}, ::Type{MOI.GreaterThan{Float64}}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{MOI.ScalarQuadraticFunction{Float64}}, ::Type{MOI.EqualTo{Float64}}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{MOI.ScalarQuadraticFunction{Float64}}, ::Type{MOI.Interval{Float64}}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{<:SF}, ::Type{<:SS}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{VAF}, ::Type{<:VLS}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{VOV}, ::Type{<:VLS}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{VAF}, ::Type{MOI.SecondOrderCone}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{VOV}, ::Type{MOI.SecondOrderCone}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{VOV}, ::Type{MOI.Complements}) = true
MOI.supports_constraint(::OnePhaseSolver, ::Type{VAF}, ::Type{MOI.Complements}) = true


########################################################
## BEGIN ModelReader CODE (with minor edits)
########################################################

#mutable struct OnePhaseMathProgModel <: MathProgBase.AbstractMathProgModel
#mutable struct OnePhaseMathProgModel
mutable struct OnePhaseMathProgModel <: MathProgBase.AbstractMathProgModel
  options
  inner::OnePhaseProblem

  eval :: Union{MOI.AbstractNLPEvaluator, Nothing}
  numVar :: Int
  numConstr :: Int
  x :: Vector{Float64}
  y :: Vector{Float64}
  lvar :: Vector{Float64}
  uvar :: Vector{Float64}
  lcon :: Vector{Float64}
  ucon :: Vector{Float64}
  sense :: Symbol
  status :: Symbol
end


#=
MOI.NonlinearModel(solver :: OnePhaseSolver) = OnePhaseMathProgModel(solver.options,OnePhaseProblem(),nothing,
                                                                   0,
                                                                   0,
                                                                   Float64[],
                                                                   Float64[],
                                                                   Float64[],
                                                                   Float64[],
                                                                   Float64[],
                                                                   Float64[],
                                                                   :Min,
                                                                   :Uninitialized)
=#
MathProgBase.NonlinearModel(solver :: OnePhaseSolver) = OnePhaseMathProgModel(solver.options,OnePhaseProblem(),nothing,
                                                                   0,
                                                                   0,
                                                                   Float64[],
                                                                   Float64[],
                                                                   Float64[],
                                                                   Float64[],
                                                                   Float64[],
                                                                   Float64[],
                                                                   :Min,
                                                                   :Uninitialized)

#println("########################################", modelTrial)
#=
MOI.RawSolver() = OnePhaseMathProgModel(this.options,OnePhaseProblem(),nothing,
                                                                   0,
                                                                   0,
                                                                   Float64[],
                                                                   Float64[],
                                                                   Float64[],
                                                                   Float64[],
                                                                   Float64[],
                                                                   Float64[],
                                                                   :Min,
                                                                   :Uninitialized)

=#
function MPB.loadproblem!(m :: OnePhaseMathProgModel,
                                   numVar, numConstr,
                                   l, u, lb, ub,
                                   sense,
                                   eval :: MOI.AbstractNLPEvaluator)

  # TODO: :JacVec is not yet available.
  # [:Grad, :Jac, :JacVec, :Hess, :HessVec, :ExprGraph]
  MPB.initialize(eval, [:Grad, :Jac, :Hess, :HessVec, :ExprGraph])
  m.numVar = numVar
  m.numConstr = numConstr
  m.x = zeros(numVar)
  m.y = zeros(numConstr)
  m.eval = eval
  m.lvar = l
  m.uvar = u
  m.lcon = lb
  m.ucon = ub
  m.sense = sense
end

MPB.setwarmstart!(m :: OnePhaseSolver, x) = (m.x = x)
MPB.status(m :: OnePhaseSolver) = m.inner.status
MPB.getsolution(m :: OnePhaseSolver) = m.inner.x
MPB.getobjval(m :: OnePhaseMathProgModel) = m.inner.iter.cache.fval
#MPB.eval_f(m.eval, m.x)

mutable struct MathProgNLPModel <: AbstractNLPModel
  meta :: NLPModelMeta
  mpmodel :: OnePhaseMathProgModel
  counters :: Counters      # Evaluation counters.

  jrows :: Vector{Int}      # Jacobian sparsity pattern.
  jcols :: Vector{Int}
  jvals :: Vector{Float64}  # Room for the constraints Jacobian.

  hrows :: Vector{Int}      # Hessian sparsity pattern.
  hcols :: Vector{Int}
  hvals :: Vector{Float64}  # Room for the Lagrangian Hessian.
end

#Utilities Begin Here
"""
    parser_SAF(fun, set, linrows, lincols, linvals, nlin, lin_lcon, lin_ucon)
Parse a `ScalarAffineFunction` fun with its associated set.
`linrows`, `lincols`, `linvals`, `lin_lcon` and `lin_ucon` are updated.
"""
function parser_SAF(fun, set, linrows, lincols, linvals, nlin, lin_lcon, lin_ucon)

  # Parse a ScalarAffineTerm{Float64}(coefficient, variable_index)
  for term in fun.terms
    push!(linrows, nlin + 1)
    push!(lincols, term.variable_index.value)
    push!(linvals, term.coefficient)
  end

  if typeof(set) in (MOI.Interval{Float64}, MOI.GreaterThan{Float64})
    push!(lin_lcon, -fun.constant + set.lower)
  elseif typeof(set) == MOI.EqualTo{Float64}
    push!(lin_lcon, -fun.constant + set.value)
  else
    push!(lin_lcon, -Inf)
  end

  if typeof(set) in (MOI.Interval{Float64}, MOI.LessThan{Float64})
    push!(lin_ucon, -fun.constant + set.upper)
  elseif typeof(set) == MOI.EqualTo{Float64}
    push!(lin_ucon, -fun.constant + set.value)
  else
    push!(lin_ucon, Inf)
  end
end

"""
    parser_VAF(fun, set, linrows, lincols, linvals, nlin, lin_lcon, lin_ucon)
Parse a `VectorAffineFunction` fun with its associated set.
`linrows`, `lincols`, `linvals`, `lin_lcon` and `lin_ucon` are updated.
"""
function parser_VAF(fun, set, linrows, lincols, linvals, nlin, lin_lcon, lin_ucon)

  # Parse a VectorAffineTerm{Float64}(output_index, scalar_term)
  for term in fun.terms
    push!(linrows, nlin + term.output_index)
    push!(lincols, term.scalar_term.variable_index.value)
    push!(linvals, term.scalar_term.coefficient)
  end

  if typeof(set) in (MOI.Nonnegatives, MOI.Zeros)
    append!(lin_lcon, -fun.constants)
  else
    append!(lin_lcon, -Inf * ones(set.dimension))
  end

  if typeof(set) in (MOI.Nonpositives, MOI.Zeros)
    append!(lin_ucon, -fun.constants)
  else
    append!(lin_ucon, Inf * ones(set.dimension))
  end
end

"""
    parser_MOI(moimodel)
Parse linear constraints of a `MOI.ModelLike`.
"""
function parser_MOI(moimodel)

  # Variables associated to linear constraints
  nlin = 0
  linrows = Int[]
  lincols = Int[]
  linvals = Float64[]
  lin_lcon = Float64[]
  lin_ucon = Float64[]

  contypes = MOI.get(moimodel, MOI.ListOfConstraints())
  for (F, S) in contypes
    F == MOI.SingleVariable && continue
    F <: AF || @warn("Function $F is not supported.")
    S <: LS || @warn("Set $S is not supported.")

    conindices = MOI.get(moimodel, MOI.ListOfConstraintIndices{F, S}())
    for cidx in conindices
      fun = MOI.get(moimodel, MOI.ConstraintFunction(), cidx)
      set = MOI.get(moimodel, MOI.ConstraintSet(), cidx)
      if typeof(fun) <: SAF
        parser_SAF(fun, set, linrows, lincols, linvals, nlin, lin_lcon, lin_ucon)
        nlin += 1
      end
      if typeof(fun) <: VAF
        parser_VAF(fun, set, linrows, lincols, linvals, nlin, lin_lcon, lin_ucon)
        nlin += set.dimension
      end
    end
  end
  coo = COO(linrows, lincols, linvals)
  nnzj = length(linvals)
  lincon = LinearConstraints(coo, nnzj)

  return nlin, lincon, lin_lcon, lin_ucon
end

"""
    parser_JuMP(jmodel)
Parse variables informations of a `JuMP.Model`.
"""
function parser_JuMP(jmodel)

  # Number of variables and bounds constraints
  nvar = Int(num_variables(jmodel))
  vars = all_variables(jmodel)
  lvar = map(
    var -> is_fixed(var) ? fix_value(var) : (has_lower_bound(var) ? lower_bound(var) : -Inf),
    vars,
  )
  uvar = map(
    var -> is_fixed(var) ? fix_value(var) : (has_upper_bound(var) ? upper_bound(var) : Inf),
    vars,
  )

  # Initial solution
  x0 = zeros(nvar)
  for (i, val) in enumerate(start_value.(vars))
    if val !== nothing
      x0[i] = val
    end
  end

  return nvar, lvar, uvar, x0
end

"""
    parser_objective_MOI(moimodel, nvar)
Parse linear and quadratic objective of a `MOI.ModelLike`.
"""
function parser_objective_MOI(moimodel, nvar)

  # Variables associated to linear and quadratic objective
  type = "UNKNOWN"
  constant = 0.0
  vect = spzeros(Float64, nvar)
  rows = Int[]
  cols = Int[]
  vals = Float64[]

  fobj = MOI.get(moimodel, MOI.ObjectiveFunction{OBJ}())

  # Single Variable
  if typeof(fobj) == SV
    type = "LINEAR"
    vect[fobj.variable.value] = 1.0
  end

  # Linear objective
  if typeof(fobj) == SAF
    type = "LINEAR"
    constant = fobj.constant
    for term in fobj.terms
      vect[term.variable_index.value] = term.coefficient
    end
  end

  # Quadratic objective
  if typeof(fobj) == SQF
    type = "QUADRATIC"
    constant = fobj.constant
    for term in fobj.affine_terms
      vect[term.variable_index.value] = term.coefficient
    end
    for term in fobj.quadratic_terms
      i = term.variable_index_1.value
      j = term.variable_index_2.value
      if i >= j
        push!(rows, i)
        push!(cols, j)
      else
        push!(cols, j)
        push!(rows, i)
      end
      push!(vals, term.coefficient)
    end
  end
  return Objective(type, constant, vect, COO(rows, cols, vals), length(vals))
end

"""
    parser_linear_expression(cmodel, nvar, F)
Parse linear expressions of type `GenericAffExpr{Float64,VariableRef}`.
"""
function parser_linear_expression(cmodel, nvar, F)

  # Variables associated to linear expressions
  rows = Int[]
  cols = Int[]
  vals = Float64[]
  constants = Float64[]

  # Linear least squares model
  nlinequ = 0
  F_is_array_of_containers = F isa Array{<:AbstractArray}
  if F_is_array_of_containers
    @objective(
      cmodel,
      Min,
      0.0 +
      0.5 *
      sum(sum(Fi^2 for Fi in FF if typeof(Fi) == GenericAffExpr{Float64, VariableRef}) for FF in F)
    )
    for FF in F, expr in FF
      if typeof(expr) == GenericAffExpr{Float64, VariableRef}
        nlinequ += 1
        for (i, key) in enumerate(expr.terms.keys)
          push!(rows, nlinequ)
          push!(cols, key.index.value)
          push!(vals, expr.terms.vals[i])
        end
        push!(constants, expr.constant)
      end
    end
  else
    @objective(
      cmodel,
      Min,
      0.0 + 0.5 * sum(Fi^2 for Fi in F if typeof(Fi) == GenericAffExpr{Float64, VariableRef})
    )
    for expr in F
      if typeof(expr) == GenericAffExpr{Float64, VariableRef}
        nlinequ += 1
        for (i, key) in enumerate(expr.terms.keys)
          push!(rows, nlinequ)
          push!(cols, key.index.value)
          push!(vals, expr.terms.vals[i])
        end
        push!(constants, expr.constant)
      end
    end
  end
  moimodel = backend(cmodel)
  lls = parser_objective_MOI(moimodel, nvar)
  return lls, LinearEquations(COO(rows, cols, vals), constants, length(vals)), nlinequ
end

"""
    parser_nonlinear_expression(cmodel, nvar, F)
Parse nonlinear expressions of type `NonlinearExpression`.
"""
function parser_nonlinear_expression(cmodel, nvar, F)

  # Nonlinear least squares model
  nnlnequ = 0
  F_is_array_of_containers = F isa Array{<:AbstractArray}
  if F_is_array_of_containers
    nnlnequ = sum(sum(typeof(Fi) == NonlinearExpression for Fi in FF) for FF in F)
    if nnlnequ > 0
      @NLobjective(
        cmodel,
        Min,
        0.5 * sum(sum(Fi^2 for Fi in FF if typeof(Fi) == NonlinearExpression) for FF in F)
      )
    end
  else
    nnlnequ = sum(typeof(Fi) == NonlinearExpression for Fi in F)
    if nnlnequ > 0
      @NLobjective(cmodel, Min, 0.5 * sum(Fi^2 for Fi in F if typeof(Fi) == NonlinearExpression))
    end
  end
  ceval = cmodel.nlp_data == nothing ? nothing : NLPEvaluator(cmodel)
  (ceval != nothing) && (nnlnequ == 0) && MOI.initialize(ceval, [:Grad, :Jac, :Hess, :HessVec])  # Add :JacVec when available
  (ceval != nothing) &&
    (nnlnequ > 0) &&
    MOI.initialize(ceval, [:Grad, :Jac, :Hess, :HessVec, :ExprGraph])  # Add :JacVec when available

  if nnlnequ == 0
    Feval = nothing
  else
    Fmodel = JuMP.Model()
    @variable(Fmodel, x[1:nvar])
    JuMP._init_NLP(Fmodel)
    @objective(Fmodel, Min, 0.0)
    Fmodel.nlp_data.user_operators = cmodel.nlp_data.user_operators
    if F_is_array_of_containers
      for FF in F, Fi in FF
        if typeof(Fi) == NonlinearExpression
          expr = ceval.subexpressions_as_julia_expressions[Fi.index]
          replace!(expr, x)
          expr = :($expr == 0)
          JuMP.add_NL_constraint(Fmodel, expr)
        end
      end
    else
      for Fi in F
        if typeof(Fi) == NonlinearExpression
          expr = ceval.subexpressions_as_julia_expressions[Fi.index]
          replace!(expr, x)
          expr = :($expr == 0)
          JuMP.add_NL_constraint(Fmodel, expr)
        end
      end
    end
    Feval = NLPEvaluator(Fmodel)
    MOI.initialize(Feval, [:Grad, :Jac, :Hess, :HessVec])  # Add :JacVec when available
  end
  return ceval, Feval, nnlnequ
end

mutable struct COO
  rows::Vector{Int}
  cols::Vector{Int}
  vals::Vector{Float64}
end

COO() = COO(Int[], Int[], Float64[])

mutable struct LinearConstraints
  jacobian::COO
  nnzj::Int
end

mutable struct LinearEquations
  jacobian::COO
  constants::Vector{Float64}
  nnzj::Int
end

mutable struct Objective
  type::String
  constant::Float64
  gradient::SparseVector{Float64}
  hessian::COO
  nnzh::Int
end

#Utilities End Here

mutable struct MathOptNLPModel <: AbstractNLPModel
  meta::NLPModelMeta
  eval::Union{MOI.AbstractNLPEvaluator, Nothing}
  lincon::LinearConstraints
  obj::Objective
  counters::Counters
end

"Construct a `MathProgNLPModel` from a `OnePhaseMathProgModel`."
function MathProgNLPModel(mpmodel :: OnePhaseMathProgModel; name :: String="Generic")

  nvar = mpmodel.numVar
  lvar = mpmodel.lvar
  uvar = mpmodel.uvar

  nlin = length(mpmodel.eval.m.linconstr)         # Number of linear constraints.
  nquad = length(mpmodel.eval.m.quadconstr)       # Number of quadratic constraints.
  nnln = length(mpmodel.eval.m.nlpdata.nlconstr)  # Number of nonlinear constraints.
  ncon = mpmodel.numConstr                        # Total number of constraints.
  lcon = mpmodel.lcon
  ucon = mpmodel.ucon

  jrows, jcols = MathProgBase.jac_structure(mpmodel.eval)
  hrows, hcols = MathProgBase.hesslag_structure(mpmodel.eval)
  nnzj = length(jrows)
  nnzh = length(hrows)

  meta = NLPModelMeta(nvar,
                      x0=mpmodel.x,
                      lvar=lvar,
                      uvar=uvar,
                      ncon=ncon,
                      y0=zeros(ncon),
                      lcon=lcon,
                      ucon=ucon,
                      nnzj=nnzj,
                      nnzh=nnzh,
                      lin=collect(1:nlin),  # linear constraints appear first in MPB
                      nln=collect(nlin+1:ncon),
                      minimize=(mpmodel.sense == :Min),
                      islp=MOI.isobjlinear(mpmodel.eval) & (nlin == ncon),
                      name=name,
                      )

  return MathProgNLPModel(meta,
                      mpmodel,
                      Counters(),
                      jrows,
                      jcols,
                      zeros(nnzj),  # jvals
                      hrows,
                      hcols,
                      zeros(nnzh),  # hvals
                      )
end

"""
    MathOptNLPModel(model, name="Generic")
Construct a `MathOptNLPModel` from a `JuMP` model.
"""
function MathOptNLPModel(jmodel::JuMP.Model; name::String = "Generic")
  ##println("HHEEEEEEEEEEEEEERRRRRRRRRRRRRRRRRREEEEEEEEEEEEEEEEEE", name)
  nvar, lvar, uvar, x0 = parser_JuMP(jmodel)

  nnln = num_nl_constraints(jmodel)

  nl_lcon = nnln == 0 ? Float64[] : map(nl_con -> nl_con.lb, jmodel.nlp_data.nlconstr)
  nl_ucon = nnln == 0 ? Float64[] : map(nl_con -> nl_con.ub, jmodel.nlp_data.nlconstr)

  eval = jmodel.nlp_data == nothing ? nothing : NLPEvaluator(jmodel)
  (eval != nothing) && MOI.initialize(eval, [:Grad, :Jac, :Hess, :HessVec])  # Add :JacVec when available

  nl_nnzj = nnln == 0 ? 0 : sum(length(nl_con.grad_sparsity) for nl_con in eval.constraints)
  nl_nnzh =
    (((eval != nothing) && eval.has_nlobj) ? length(eval.objective.hess_I) : 0) +
    (nnln == 0 ? 0 : sum(length(nl_con.hess_I) for nl_con in eval.constraints))

  moimodel = backend(jmodel)
  nlin, lincon, lin_lcon, lin_ucon = parser_MOI(moimodel)

  if (eval != nothing) && eval.has_nlobj
    obj = Objective("NONLINEAR", 0.0, spzeros(Float64, nvar), COO(), 0)
  else
    obj = parser_objective_MOI(moimodel, nvar)
  end

  ncon = nlin + nnln
  lcon = vcat(lin_lcon, nl_lcon)
  ucon = vcat(lin_ucon, nl_ucon)
  nnzj = lincon.nnzj + nl_nnzj
  nnzh = obj.nnzh + nl_nnzh

  meta = NLPModelMeta(
    nvar,
    x0 = x0,
    lvar = lvar,
    uvar = uvar,
    ncon = ncon,
    nlin = nlin,
    nnln = nnln,
    y0 = zeros(ncon),
    lcon = lcon,
    ucon = ucon,
    nnzj = nnzj,
    nnzh = nnzh,
    lin = collect(1:nlin),
    nln = collect((nlin + 1):ncon),
    minimize = objective_sense(jmodel) == MOI.MIN_SENSE,
    islp = (obj.type == "LINEAR") && (nnln == 0),
    name = name,
  )

  return MathOptNLPModel(meta, eval, lincon, obj, Counters())
end

import Base.show
show(nlp :: MathProgNLPModel) = show(nlp.mpmodel)

function obj(nlp :: MathProgNLPModel, x :: Array{Float64})
  NLPModels.increment!(nlp, :neval_obj)
  return MathProgBase.eval_f(nlp.mpmodel.eval, x)
end

function obj(nlp :: MathOptNLPModel, x :: Array{Float64})
  NLPModels.increment!(nlp, :neval_obj)
  return MOI.eval_objective(nlp.eval, x)
end

function grad(nlp :: MathProgNLPModel, x :: Array{Float64})
  g = zeros(nlp.meta.nvar)
  return grad!(nlp, x, g)
end

function grad!(nlp :: MathProgNLPModel, x :: Array{Float64}, g :: Array{Float64})
  NLPModels.increment!(nlp, :neval_grad)
  MathProgBase.eval_grad_f(nlp.mpmodel.eval, g, x)
  return g
end

function grad(nlp :: MathOptNLPModel, x :: Array{Float64})
  g = zeros(nlp.meta.nvar)
  return grad!(nlp, x, g)
end

function grad!(nlp :: MathOptNLPModel, x :: Array{Float64}, g :: Array{Float64})
  NLPModels.increment!(nlp, :neval_grad)
  MOI.eval_objective_gradient(nlp.eval, g, x)
  return g
end

function cons(nlp :: MathProgNLPModel, x :: Array{Float64})
  c = zeros(nlp.meta.ncon)
  return cons!(nlp, x, c)
end

function cons(nlp :: MathOptNLPModel, x :: Array{Float64})
  c = zeros(nlp.meta.ncon)
  return cons!(nlp, x, c)
end

function cons!(nlp :: MathProgNLPModel, x :: Array{Float64}, c :: Array{Float64})
  NLPModels.increment!(nlp, :neval_cons)
  MathProgBase.eval_g(nlp.mpmodel.eval, c, x)
  return c
end

function cons!(nlp :: MathOptNLPModel, x :: Array{Float64}, c :: Array{Float64})
  NLPModels.increment!(nlp, :neval_cons)
  MOI.eval_constraint(nlp.eval, c, x)
  return c
end


function jac_coord(nlp :: MathProgNLPModel, x :: Array{Float64})
  NLPModels.increment!(nlp, :neval_jac)
  MOI.eval_jac_g(nlp.mpmodel.eval, nlp.jvals, x)
  return (nlp.jrows, nlp.jcols, nlp.jvals)
end

function jac_coord(nlp :: MathOptNLPModel, x :: Array{Float64})
  NLPModels.increment!(nlp, :neval_jac)
  MOI.eval_constraint_jacobian(nlp.eval, nlp.lincon.jacobian.vals, x)
  return (nlp.lincon.jacobian.rows, nlp.lincon.jacobian.cols, nlp.lincon.jacobian.vals)
  #=
  increment!(nlp, :neval_jac)
  if nlp.meta.nlin > 0
    vals[1:(nlp.lincon.nnzj)] .= nlp.lincon.jacobian.vals[1:(nlp.lincon.nnzj)]
  end
  if nlp.meta.nnln > 0
    MOI.eval_constraint_jacobian(nlp.eval, view(vals, (nlp.lincon.nnzj + 1):(nlp.meta.nnzj)), x)
  end
  return vals
  =#
end

function NLPModels.jac_coord!(nlp::MathOptNLPModel, x::AbstractVector, vals::AbstractVector)
  increment!(nlp, :neval_jac)
  if nlp.meta.nlin > 0
    vals[1:(nlp.lincon.nnzj)] .= nlp.lincon.jacobian.vals[1:(nlp.lincon.nnzj)]
  end
  if nlp.meta.nnln > 0
    MOI.eval_constraint_jacobian(nlp.eval, view(vals, (nlp.lincon.nnzj + 1):(nlp.meta.nnzj)), x)
  end
  return vals
end

function jac(nlp :: MathOptNLPModel, x :: Array{Float64})
  return SparseArrays.sparse(jac_coord(nlp, x)..., nlp.meta.ncon, nlp.meta.nvar)
end

function jac(nlp :: MathProgNLPModel, x :: Array{Float64})
  return SparseArrays.sparse(jac_coord(nlp, x)..., nlp.meta.ncon, nlp.meta.nvar)
end


function jprod(nlp :: MathProgNLPModel, x :: Array{Float64}, v :: Array{Float64})
  Jv = zeros(nlp.meta.ncon)
  return jprod!(nlp, x, v, Jv)
end

function jprod!(nlp :: MathProgNLPModel,
                x :: Array{Float64},
                v :: Array{Float64},
                Jv :: Array{Float64})
  nlp.counters.neval_jac -= 1
  NLPModels.increment!(nlp, :neval_jprod)
  Jv[:] = jac(nlp, x) * v
  return Jv
end

function jtprod(nlp :: MathProgNLPModel, x :: Array{Float64}, v :: Array{Float64})
  Jtv = zeros(nlp.meta.nvar)
  return jtprod!(nlp, x, v, Jtv)
end

function jtprod!(nlp :: MathProgNLPModel,
                x :: Array{Float64},
                v :: Array{Float64},
                Jtv :: Array{Float64})
  nlp.counters.neval_jac -= 1
  NLPModels.increment!(nlp, :neval_jtprod)
  Jtv[1:nlp.meta.nvar] = jac(nlp, x)' * v
  return Jtv
end

# Uncomment if/when :JacVec becomes available in MPB.
# "Evaluate the Jacobian-vector product at `x`."
# function jprod(nlp :: MathProgNLPModel, x :: Array{Float64}, v :: Array{Float64})
#   jv = zeros(nlp.meta.ncon)
#   return jprod!(nlp, x, v, jv)
# end
#
# "Evaluate the Jacobian-vector product at `x` in place."
# function jprod!(nlp :: MathProgNLPModel, x :: Array{Float64}, v :: Array{Float64}, jv :: Array{Float64})
#   NLPModels.increment!(nlp, :neval_jprod)
#   MOI.eval_jac_prod(nlp.mpmodel.eval, jv, x, v)
#   return jv
# end
#
# "Evaluate the transposed-Jacobian-vector product at `x`."
# function jtprod(nlp :: MathProgNLPModel, x :: Array{Float64}, v :: Array{Float64})
#   jtv = zeros(nlp.meta.nvar)
#   return jtprod!(nlp, x, v, jtv)
# end
#
# "Evaluate the transposed-Jacobian-vector product at `x` in place."
# function jtprod!(nlp :: MathProgNLPModel, x :: Array{Float64}, v :: Array{Float64}, jtv :: Array{Float64})
#   NLPModels.increment!(nlp, :neval_jtprod)
#   MathProgBase.eval_jac_prod_t(nlp.mpmodel.eval, jtv, x, v)
#   return jtv
# end

function hess_coord(nlp :: MathProgNLPModel, x :: Array{Float64};
    obj_weight :: Float64=1.0, y :: Array{Float64}=zeros(nlp.meta.ncon))
  NLPModels.increment!(nlp, :neval_hess)
  MOI.eval_hesslag(nlp.mpmodel.eval, nlp.hvals, x, obj_weight, y)
  return (nlp.hrows, nlp.hcols, nlp.hvals)
end

function hess(nlp :: MathProgNLPModel, x :: Array{Float64};
    obj_weight :: Float64=1.0, y :: Array{Float64}=zeros(nlp.meta.ncon))
  return SparseArrays.sparse(hess_coord(nlp, x, y=y, obj_weight=obj_weight)..., nlp.meta.nvar, nlp.meta.nvar)
end

function hess_coord(nlp :: MathOptNLPModel, x :: Array{Float64};
    obj_weight :: Float64=1.0, y :: Array{Float64}=zeros(nlp.meta.ncon))
  NLPModels.increment!(nlp, :neval_hess)
  MOI.eval_hessian_lagrangian(nlp.eval, nlp.obj.hessian.vals, x, obj_weight, y)
  return (nlp.obj.hessian.rows, nlp.obj.hessian.cols, nlp.obj.hessian.vals)
end

function hess(nlp :: MathOptNLPModel, x :: Array{Float64};
    obj_weight :: Float64=1.0, y :: Array{Float64}=zeros(nlp.meta.ncon))
  return SparseArrays.sparse(hess_coord(nlp, x, y=y, obj_weight=obj_weight)..., nlp.meta.nvar, nlp.meta.nvar)
end

function hprod(nlp :: MathProgNLPModel, x :: Array{Float64}, v :: Array{Float64};
    obj_weight :: Float64=1.0, y :: Array{Float64}=zeros(nlp.meta.ncon))
  hv = zeros(nlp.meta.nvar)
  return hprod!(nlp, x, v, hv, obj_weight=obj_weight, y=y)
end

function hprod!(nlp :: MathProgNLPModel, x :: Array{Float64}, v :: Array{Float64},
    hv :: Array{Float64};
    obj_weight :: Float64=1.0, y :: Array{Float64}=zeros(nlp.meta.ncon))
  NLPModels.increment!(nlp, :neval_hprod)
  MOI.eval_hesslag_prod(nlp.mpmodel.eval, hv, x, v, obj_weight, y)
  return hv
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

function MOI.optimize!(m :: OnePhaseMathProgModel)
    t = time()
    nlp = MathProgNLPModel(m)

    pars = create_pars_JuMP(m.options)

    iter, status, hist, t, err, timer = one_phase_solve(nlp,pars)

    m.inner.status = status_One_Phase_To_JuMP(status)
    m.inner.x = get_original_x(iter)
    m.inner.obj_val = iter.cache.fval
    m.inner.lambda = get_y(iter)
    m.inner.solve_time = time() - t

    # custom one phase features
    m.inner.pars = pars
    m.inner.iter = iter
    m.inner.hist = hist
end
#=
MOI.getconstrsolution(m::OnePhaseMathProgModel) = m.inner.g
MOI.getrawsolver(m::OnePhaseMathProgModel) = m.inner
MOI.getsolvetime(m::OnePhaseMathProgModel) = m.inner.solve_time

function MOI.getreducedcosts(m::OnePhaseMathProgModel)
    return get_reducedcosts(m.inner.iter)
end

function MOI.getconstrduals(m::OnePhaseMathProgModel)
    return get_constrduals(m.inner.iter)
end
=#
#setvartype!(m::OnePhaseMathProgModel, typ::Vector{Symbol}) =
#    (m.varType = map(t->rev_var_type_map[t], typ))
#=
function MOI.freemodel!(m::OnePhaseMathProgModel)
    # TO DO
end
=#

function free(model::OnePhaseMathProgModel)
    if model.inner != nothing
        KN_free(model.inner)
    end
end

function MOI.empty!(model::OnePhaseSolver)
end