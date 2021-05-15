using SparseArrays
@compat abstract type abstract_nlp end

mutable struct Class_cache
    fval::Float64
    cons::Array{Float64,1}
    grad::Array{Float64,1}
    J::SparseMatrixCSC{Float64,Int64}
    J_T::SparseMatrixCSC{Float64,Int64}
    H::SparseMatrixCSC{Float64,Int64}
    fval_updated::Bool
    cons_updated::Bool
    grad_updated::Bool
    J_updated::Bool
    H_updated::Bool

    function Class_cache()
      return new()
    end
end

function is_updated(c::Class_cache)
    return c.fval_updated && c.cons_updated && c.grad_updated && c.J_updated && c.H_updated
end

function is_updated_correction(c::Class_cache)
    return c.fval_updated && c.cons_updated && c.grad_updated && c.J_updated
end

mutable struct Class_local_info
    delta::Float64
    radius::Float64
    prev_step::Symbol
    function Class_local_info()
      return new(0.0,0.0,:blah)
    end
end


mutable struct Class_iterate
    primal_residual_intial::Array{Float64,1}
    point::Class_point
    nlp::abstract_nlp
    cache::Class_cache
    local_info::Class_local_info
    ncon::Int64
    nvar::Int64
    frac_bd::Array{Float64,1}
    frac_bd_predict::Array{Float64,1}

    # regularizer info
    a_norm_penalty_par::Float64

    function Class_iterate()
        return new();
    end

    function Class_iterate(intial_point::Class_point, nlp::abstract_nlp, local_info::Class_local_info, timer::class_advanced_timer, pars::Class_parameters)
      this = new() # zeros(length(intial_point.x)), intial_point, nlp, Class_cache(), local_info)
      this.primal_residual_intial = zeros(length(intial_point.x))
      this.point = intial_point
      this.nlp = nlp
      this.local_info = local_info
      this.nvar = length(intial_point.x)
      this.ncon = length(intial_point.y)
      this.a_norm_penalty_par = pars.a_norm_penalty;
      this.frac_bd = pars.ls.fraction_to_boundary * ones(this.ncon)
      this.frac_bd_predict = pars.ls.fraction_to_boundary_predict * ones(this.ncon)

      @assert(this.ncon == length(intial_point.y))

      this.cache = Class_cache()
      update!(this, timer, pars)
      #update_cons!(this, timer, pars)
      #update_obj!(this, timer, pars)
      #update_H!(this, timer, pars)
      #update_J!(this, timer, pars)

      #update_primal_cache!(this, timer)

      this.primal_residual_intial = get_primal_res(this);
      return this
    end
end

function get_original_x(it::Class_iterate)
     _cute_x(it.nlp, it.point.x)
end

function get_constrduals(it::Class_iterate)
    return get_constrduals(it.nlp,it.point.y)
end

function get_reducedcosts(it::Class_iterate)
    return get_reducedcosts(it.nlp,it.point.y)
end

function validate(it::Class_iterate)
    @assert(it.nvar == length(it.point.x))
    @assert(it.ncon == length(it.point.y))
    @assert(it.ncon == length(it.point.s))
end

import Base.copy

function copy(it::Class_iterate, timer::class_advanced_timer)
   validate(it)

   new_it = Class_iterate()
   new_it.point = deepcopy(it.point)

   new_it.primal_residual_intial = it.primal_residual_intial
   new_it.nlp = it.nlp
   new_it.local_info = it.local_info
   new_it.nvar = it.nvar
   new_it.ncon = it.ncon

   new_it.a_norm_penalty_par = it.a_norm_penalty_par;
   new_it.frac_bd = it.frac_bd
   new_it.frac_bd_predict = it.frac_bd_predict

   new_it.cache = it.cache
   copy_cache!(new_it, timer)

   return new_it
end

function finite_diff_check(it::Class_iterate)
   approx_grad = Calculus.gradient(it.nlp.eval_f, it.point.x)
   actual_grad = it.nlp.eval_grad_lag(it.point.x, 0.0, zeros(length(it.point.y)))
   @assert(norm(approx_grad -  actual_grad, 2) < 1e-3)

   actual_jac = it.nlp.eval_jac(it.point.x)
   for i = 1:length(it.point.s)
       a_i = x::Vector -> it.nlp.eval_a(x)[i]
       approx_grad_a_i = Calculus.gradient(a_i, it.point.x)
       @assert(norm(actual_jac[i,:] - approx_grad_a_i,2) < 1e-3)
   end
end

function feasible_ls(it::Class_iterate, dir::Class_point, pars::Class_parameters, timer::class_advanced_timer)

ccifg(io_err, n, icon, x, ci, gci, grad)

end


function dim(it::Class_iterate)
    return length(it.point.x)
end

function ncon(it::Class_iterate)
    @assert(length(length(it.point.s)) == length(length(it.point.y)))
    return length(it.point.s)
end

function get_delta(it::Class_iterate)
    return it.local_info.delta
end

function set_delta(it::Class_iterate, delta::Float64)
    it.local_info.delta = delta;
end

function get_radius(it::Class_iterate)
    return it.local_info.radius
end

function set_prev_step(it::Class_iterate, prev_step::Symbol)
    it.local_info.prev_step = prev_step;
end

function get_prev_step(it::Class_iterate)
    return it.local_info.prev_step
end
function set_radius(it::Class_iterate, radius::Float64)
    it.local_info.radius = radius;
end

function get_x(it::Class_iterate)
    return it.point.x
end

function get_y(it::Class_iterate)
    return it.point.y
end

function get_s(it::Class_iterate)
    return it.point.s
end

function get_mu(it::Class_iterate)
    return it.point.mu
end

function copy_cache!(it::Class_iterate, timer::class_advanced_timer)
    nlp = it.nlp
    start_advanced_timer(timer, "CACHE")
    start_advanced_timer(timer, "CACHE/copy")
    cache = it.cache

    new_cache = Class_cache()
    new_cache.J = deepcopy(cache.J)
    new_cache.J_T = deepcopy(cache.J_T)
    new_cache.J_updated = false
    new_cache.H = deepcopy(cache.H)
    new_cache.H_updated = false
    new_cache.cons = deepcopy(cache.cons)
    new_cache.cons_updated = false
    new_cache.fval = deepcopy(cache.fval)
    new_cache.fval_updated = false
    new_cache.grad = deepcopy(cache.grad)
    new_cache.grad_updated = false


    it.cache = new_cache

    pause_advanced_timer(timer, "CACHE/copy")
    pause_advanced_timer(timer, "CACHE")
end

function update!(it::Class_iterate, timer::class_advanced_timer, pars::Class_parameters)
    success = update_cons!(it, timer, pars) && update_obj!(it, timer, pars) && update_H!(it, timer, pars) && update_grad!(it, timer, pars) && update_J!(it, timer, pars)
    return success
end

function getbad(something::Array{Float64,1})
    for val in something
      if isbad(val)
        return val
      end
    end

    return 0.0
end

function update_cons!(it::Class_iterate, timer::class_advanced_timer, pars::Class_parameters)
    start_advanced_timer(timer, "CACHE")
    start_advanced_timer(timer, "CACHE/update_cons")
    it.cache.cons = eval_a(it.nlp, it.point.x)
    it.cache.cons_updated = true
    if isbad(it.cache.cons)
        println("NaN or Inf in cons")
        if pars.throw_error_nans
          throw(Eval_NaN_error(getbad(it.cache.cons), it.point.x,"cons"))
        else
          no_nan = false
        end
    else
      no_nan = true
    end
    pause_advanced_timer(timer, "CACHE/update_cons")
    pause_advanced_timer(timer, "CACHE")

    return no_nan
end

function update_obj!(it::Class_iterate, timer::class_advanced_timer, pars::Class_parameters)
    start_advanced_timer(timer, "CACHE")
    start_advanced_timer(timer, "CACHE/update_obj")
    it.cache.fval = eval_f(it.nlp, it.point.x)
    it.cache.fval_updated = true
    if isbad(it.cache.fval)
        println("$(it.cache.fval) in obj")
        if pars.throw_error_nans
          throw(Eval_NaN_error(it.cache.fval, it.point.x, "obj"))
        else
          no_nan = false
        end
    else
      no_nan = true
    end
    pause_advanced_timer(timer, "CACHE/update_obj")
    pause_advanced_timer(timer, "CACHE")

    return no_nan
end

function update_H!(it::Class_iterate, timer::class_advanced_timer, pars::Class_parameters)
    start_advanced_timer(timer, "CACHE")
    start_advanced_timer(timer, "CACHE/update_H")

    # regularizer
    begin

      lambda = it.point.mu #* use_prox(it)
      y_rx = lambda * it.a_norm_penalty_par

      x = it.point.x

      it.cache.H = eval_lag_hess(it.nlp, x, it.point.y + y_rx, 1.0)
      #it.cache.H = eval_lag_hess(it.nlp, it.point.x, it.point.y, 1.0)
      it.cache.H_updated = true
    end

    nzH = nonzeros(it.cache.H)
    if isbad(nzH)
        println("NaN or Inf in H")
        if pars.throw_error_nans
          throw(Eval_NaN_error(getbad(nzH), it.point.x, "H"))
        else
          no_nan = false
        end
    else
      no_nan = true
    end
    pause_advanced_timer(timer, "CACHE/update_H")
    pause_advanced_timer(timer, "CACHE")

    return no_nan
end

function update_grad!(it::Class_iterate, timer::class_advanced_timer, pars::Class_parameters)
    start_advanced_timer(timer, "CACHE")
    start_advanced_timer(timer, "CACHE/update_grad")
    it.cache.grad = eval_grad_f(it.nlp, it.point.x)
    it.cache.grad_updated = true
    if isbad(it.cache.grad)
        println("NaN or Inf in grad")
        if pars.throw_error_nans
          throw(Eval_NaN_error(getbad(it.cache.grad), it.point.x, "grad"))
        else
          no_nan = false
        end
    else
      no_nan = true
    end
    pause_advanced_timer(timer, "CACHE/update_grad")
    pause_advanced_timer(timer, "CACHE")

    return no_nan
end

function update_J!(it::Class_iterate, timer::class_advanced_timer, pars::Class_parameters)
    start_advanced_timer(timer, "CACHE")
    start_advanced_timer(timer, "CACHE/update_J")
    @assert(length(it.point.x) == length(it.cache.grad))
    it.cache.J = eval_jac(it.nlp, it.point.x)
    it.cache.J_T = it.cache.J'
    it.cache.J_updated = true
    if isbad(nonzeros(it.cache.J))
        println("NaN or Inf in J")
        if pars.throw_error_nans
          throw(Eval_NaN_error(getbad(it.cache.J), it.point.x, "J"))
        else
          no_nan = false
        end
    else
      no_nan = true
    end
    pause_advanced_timer(timer, "CACHE/update_J")
    pause_advanced_timer(timer, "CACHE")

    return no_nan
end
