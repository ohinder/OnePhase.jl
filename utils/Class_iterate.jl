abstract abstract_nlp;

type Class_cache
    fval::Float64
    cons::Array{Float64,1}
    grad::Array{Float64,1}
    J::SparseMatrixCSC{Float64,Int32}
    H::SparseMatrixCSC{Float64,Int32}
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

type Class_local_info
    delta::Float64
    radius::Float64
    prev_step::Symbol
end


type Class_iterate
    primal_residual_intial::Array{Float64,1}
    point::Class_point
    nlp::abstract_nlp
    cache::Class_cache
    local_info::Class_local_info
    ncon::Int64
    nvar::Int64

    function Class_iterate(intial_point::Class_point, nlp::abstract_nlp, local_info::Class_local_info, timer::class_advanced_timer, pars::Class_parameters)
      this = new(zeros(length(intial_point.x)), intial_point, nlp, Class_cache(), local_info)
      this.nvar = length(intial_point.x)
      this.ncon = length(intial_point.y)
      @assert(this.ncon == length(intial_point.y))

      this.cache = Class_cache()
      update_grad!(this, timer, pars)
      update_cons!(this, timer, pars)
      update_obj!(this, timer, pars)
      update_H!(this, timer, pars)
      update_J!(this, timer, pars)

      #update_primal_cache!(this, timer)

      this.primal_residual_intial = get_primal_res(this);
      return this
    end

    function Class_iterate()
      return new()
    end
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

   new_it.cache = it.cache
   copy_cache!(new_it, timer)

   return new_it
end

function finite_diff_check(it::Class_iterate)
   approx_grad = Calculus.gradient(it.nlp.eval_f, it.point.x)
   actual_grad = it.nlp.eval_grad_lag(it.point.x, zeros(length(it.point.y)))
   @assert(norm(approx_grad -  actual_grad, 2) < 1e-3)

   actual_jac = it.nlp.eval_jac(it.point.x)
   for i = 1:length(it.point.s)
       a_i = x::Vector -> it.nlp.eval_a(x)[i]
       approx_grad_a_i = Calculus.gradient(a_i, it.point.x)
       @assert(norm(actual_jac[i,:] - approx_grad_a_i,2) < 1e-3)
   end
end

function dual_bounds(it::Class_iterate, y::Array{Float64,1}, dy::Array{Float64,1}, comp_feas::Float64)
    lb = 0.0
    ub = 1.0
    mu = it.point.mu
    for i = 1:length(y)
        s = it.point.s[i]
        @assert(s > 0.0)

        safety_factor = 1.01
        safety_add = 0.0

        #ub_dyi = (mu / (comp_feas * s) - y[i]) / dy[i]
        #lb_dyi = (mu * comp_feas / s - y[i]) / dy[i]
        #@show mu * comp_feas / s, y[i]

        ub_dyi = (mu / (comp_feas * s * dy[i]) - y[i] / dy[i])
        lb_dyi = (mu * comp_feas / (s * dy[i]) - y[i]  / dy[i])

        if dy[i] > 0.0
          lb = max( lb_dyi * safety_factor + safety_add, lb)
          ub = min( ub_dyi / safety_factor - safety_add, ub)

          #if lb < ub
          #  @assert( y[i] + ub * dy[i] <= mu / (comp_feas * s) )
          #  @assert( y[i] + lb * dy[i] >= comp_feas * mu / s )
          #end
        elseif dy[i] < 0.0
          lb = max(ub_dyi * safety_factor + safety_add, lb)
          ub = min(lb_dyi / safety_factor - safety_add, ub)

          #if lb < ub
          #  @assert( y[i] + lb * dy[i] <= mu / (comp_feas * s) )
          #  @assert( y[i] + ub * dy[i] >= comp_feas * mu / s )
          #end
        elseif lb_dyi >= 0.0 || ub_dyi <= 0.0
          lb = 0.0
          ub = -1.0
        end

        if lb < ub
          mu_lb = (y[i] + lb * dy[i]) * s
          mu_ub = (y[i] + ub * dy[i]) * s

          if mu_lb > mu / comp_feas
              println("HERE1")
          end

          if mu_ub > mu / comp_feas
            println("HERE2")
          end

          if mu_ub < mu * comp_feas
            println("HERE3")
          end

          if mu_lb < mu * comp_feas
            println("HERE4")
          end
        end

        #|| mu_lb < mu * comp_feas
          #@show mu_lb, mu_ub
          #@show comp_feas * mu, mu / comp_feas
          #@show lb
          #@show mu, dy[i]
          #@assert(mu_lb <= mu / comp_feas)
          #@assert(mu_ub <= mu / comp_feas)
          #@assert(mu_lb >= comp_feas * mu)
          #@assert(mu_ub >= comp_feas * mu)
        #end
    end

    return lb, ub
end

function feasible_ls(it::Class_iterate, dir::Class_point, pars::Class_parameters, timer::class_advanced_timer)

ccifg(io_err, n, icon, x, ci, gci, grad)

end


function move(it::Class_iterate, dir::Class_point, step_size::Float64, pars::Class_parameters, timer::class_advanced_timer)
    if pars.move_type == :primal_dual
      new_it = copy(it, timer)
      new_it.point.x += dir.x * step_size
      update_cons!(new_it, timer, pars)

      new_a = get_cons(new_it)

      #nl_s = new_it.point.s + step_size * dir.s
      new_it.point.primal_scale += dir.primal_scale * step_size
      nl_s = new_a - new_it.point.primal_scale * it.primal_residual_intial



      if pars.s_update == :careful
        new_it.point.s = nl_s
      elseif pars.s_update == :loose
        rf = 1.5
        l_s = new_it.point.s + step_size * dir.s
        new_it.point.s = min(nl_s * rf, max(nl_s / rf, l_s))
      else
        error("pars.s_update option $(pars.s_update) not avaliable")
      end

      if !all(new_it.point.s .>= it.point.s * pars.fraction_to_boundary)
          return it, :s_bound, step_size
      end

      if pars.mu_update == :static || (pars.mu_update == :dynamic_agg && dir.mu == 0.0)
        new_it.point.mu += dir.mu * step_size
      elseif pars.mu_update == :dynamic || pars.mu_update == :dynamic_agg
        new_it.point.mu = mean(new_it.point.s .* (new_it.point.y + step_size * dir.y))
      else
        error("this mu_update option is not valid")
      end

      if minimum(new_it.point.s) > 0.0
        lb1, ub1 = dual_bounds(new_it, new_it.point.y, dir.y, pars.comp_feas)
        lb2, ub2 = dual_bounds(new_it, new_it.point.y, dir.y, pars.comp_feas_agg)

        if pars.move_primal_seperate_to_dual
            if lb1 < ub1
                theta = (1.0 - pars.fraction_to_boundary)
                step_size_D_boundary = theta / maximum([theta; -dir.y ./ it.point.y])
                step_size_D_max = min(step_size_D_boundary, ub1)

                if pars.dual_ls
                  scale_D = dual_scale(new_it, pars)
                  scale_mu = dual_scale(new_it, pars)
                  #scale = 1.0;
                  ∇a = get_jac(new_it)
                  q = [scale_D * ∇a' * dir.y; scale_mu * new_it.point.s .* dir.y];
                  predicted_dual_res = step_size * (hess_product(it, dir.x) + get_delta(it) * dir.x) + eval_grad_lag(it)
                  res = [scale_D * predicted_dual_res; -scale_mu * comp(new_it)]
                  step_size_D = sum(res .* q) / sum(q.^2)
                  step_size_D = max(lb1,min(step_size_D,step_size_D_max))
                  #=@show sum(q.^2)
                  @show norm(dir.x, 2)
                  @show get_delta(it),  norm(eval_grad_lag(it)), norm(new_it.point.s,2)=#
                else
                  step_size_D = step_size_D_max
                end

                new_it.point.y += dir.y * step_size_D

                #@show dir.y
                #@show step_size_D

                #if !all(new_it.point.y .>= it.point.y * pars.fraction_to_boundary)
                #    return it, false, step_size_D
                #end

                #update_H!(it, timer)

                if !is_feasible(new_it, pars.comp_feas)
                  @show minimum(new_it.point.y), maximum(new_it.point.s)
                  println("line=",@__LINE__, "file=",@__FILE__)
                  my_warn("infeasibility should have been detected earlier!!!")
                  @show step_size, step_size_D
                  @show (lb1, ub1), (lb2, ub2)
                  return it, :dual_infeasible, step_size_D
                else
                  return new_it, :success, step_size_D
                end
            else
              return it, :dual_infeasible, step_size
            end
        else
          if lb1 <= step_size && step_size <= ub1
            new_it.point.y += dir.y * step_size

            if !all(new_it.point.y .>= it.point.y * pars.fraction_to_boundary)
                return it, false, step_size
            end

            if !is_feasible(new_it, pars.comp_feas)
              return it, false, step_size
            else
              return new_it, true, step_size
            end
          else
            return it, false, step_size
          end
        end
      else
          return it, false, step_size
      end
      #new_it.point.y = new_it.point.mu ./ new_it.point.s
    elseif pars.move_type == :primal
      new_it = copy(it)
      new_it.point.mu += dir.mu * step_size
      new_it.point.x += dir.x * step_size
      new_a = eval_a(new_it)

      new_it.point.primal_scale += dir.primal_scale * step_size
      new_it.point.s = new_a - new_it.point.primal_scale * it.primal_residual_intial
      #new_it.point.s += dir.s
      new_it.point.y = new_it.point.mu ./ new_it.point.s

      if is_feasible(new_it, pars.comp_feas)
        return new_it, true, step_size
      else
        return it, false, step_size
      end
    end
end
#=
function move(it::Class_iterate, dir::Class_point, step_size::Float64, pars::Class_parameters, timer::class_advanced_timer)
    if pars.move_type == :primal_dual
      new_it = copy(it)
      new_it.point.x += dir.x * step_size
      update_primal_cache!(new_it, timer)
      new_a = get_cons(new_it)

      #nl_s = new_it.point.s + step_size * dir.s
      new_it.point.primal_scale += dir.primal_scale * step_size
      nl_s = new_a - new_it.point.primal_scale * it.primal_residual_intial

      if pars.s_update == :careful
        new_it.point.s = nl_s
      elseif pars.s_update == :loose
        rf = 1.5
        l_s = new_it.point.s + step_size * dir.s
        new_it.point.s = min(nl_s * rf, max(nl_s / rf, l_s))
      else
        error("pars.s_update option $(pars.s_update) not avaliable")
      end

      if pars.mu_update == :static || (pars.mu_update == :dynamic_agg && dir.mu == 0.0)
        new_it.point.mu += dir.mu * step_size
      elseif pars.mu_update == :dynamic || pars.mu_update == :dynamic_agg
        new_it.point.mu = mean(new_it.point.s .* (new_it.point.y + step_size * dir.y))
      else
        error("this mu_update option is not valid")
      end

      if minimum(new_it.point.s) > 0.0
        lb1, ub1 = dual_bounds(it, new_it.point.y, dir.y, pars.comp_feas)
        lb2, ub2 = dual_bounds(it, new_it.point.y, dir.y, pars.comp_feas_agg)

        if pars.move_primal_seperate_to_dual
            if lb1 < ub1 / 1.01
                ∇a = get_jac(new_it)
                scale = dual_scale(new_it, pars)

                if lb2 < ub2
                    q = [scale * ∇a' * dir.y; new_it.point.s .* dir.y];
                    res = [scale * eval_grad_lag(new_it); -comp(new_it)]
                    lb = lb2
                    ub = ub2
                else
                    q = new_it.point.s .* dir.y;
                    res = -comp(new_it);
                    lb = lb1
                    ub = ub1
                end

                step_size_D = sum(res .* q) / sum(q.^2)
                step_size_D = min(ub, max(lb, step_size_D))


                new_it.point.y += dir.y * step_size_D

                #update_H!(it, timer)

                if !is_feasible(new_it, pars.comp_feas)
                  @show minimum(new_it.point.y), maximum(new_it.point.s)
                  println("line=",@__LINE__, "file=",@__FILE__)
                  my_warn("infeasibility should have been detected earlier!!!")
                  @show (lb1, ub1), (lb2, ub2)
                  return it, false, step_size_D
                else
                  return new_it, true, step_size_D
                end
            else
              return it, false, step_size
            end
        else
          if lb1 <= step_size && step_size <= ub1
            new_it.point.y += dir.y * step_size

            return new_it, true, step_size
          else
            return it, false, step_size
          end
        end
      else
          return it, false, step_size
      end
      #new_it.point.y = new_it.point.mu ./ new_it.point.s
    elseif pars.move_type == :primal
      new_it = copy(it)
      new_it.point.mu += dir.mu * step_size
      new_it.point.x += dir.x * step_size
      new_a = eval_a(new_it)

      new_it.point.primal_scale += dir.primal_scale * step_size
      new_it.point.s = new_a - new_it.point.primal_scale * it.primal_residual_intial
      #new_it.point.s += dir.s
      new_it.point.y = new_it.point.mu ./ new_it.point.s

      if is_feasible(new_it, pars.comp_feas)
        return new_it, true, step_size
      else
        return it, false, step_size
      end
    end
end=#

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
    start_advanced_timer(timer, "CACHE/copy")
    cache = it.cache

    new_cache = Class_cache()
    new_cache.J = copy(cache.J)
    new_cache.J_updated = false
    new_cache.H = copy(cache.H)
    new_cache.H_updated = false
    new_cache.cons = copy(cache.cons)
    new_cache.cons_updated = false
    new_cache.fval = copy(cache.fval)
    new_cache.fval_updated = false
    new_cache.grad = copy(cache.grad)
    new_cache.grad_updated = false


    it.cache = new_cache

    pause_advanced_timer(timer, "CACHE/copy")
end

function update!(it::Class_iterate, timer::class_advanced_timer, pars::Class_parameters)
    update_cons!(it, timer, pars)
    update_obj!(it, timer, pars)
    update_H!(it, timer, pars)
    update_grad!(it, timer, pars)
    update_J!(it, timer, pars)
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
    start_advanced_timer(timer, "CACHE/update_cons")
    it.cache.cons = eval_a(it.nlp, it.point.x)
    it.cache.cons_updated = true
    if isbad(it.cache.cons)
        println("NaN or Inf")
        throw(Eval_NaN_error(getbad(it.cache.cons), it.point.x,"cons"))
    end
    pause_advanced_timer(timer, "CACHE/update_cons")
end

function update_obj!(it::Class_iterate, timer::class_advanced_timer, pars::Class_parameters)
    start_advanced_timer(timer, "CACHE/update_obj")
    it.cache.fval = eval_f(it.nlp, it.point.x)
    it.cache.fval_updated = true
    if isbad(it.cache.fval)
        println("NaN or Inf")
        throw(Eval_NaN_error(it.cache.fval, it.point.x, "obj"))
    end
    pause_advanced_timer(timer, "CACHE/update_obj")
end

function update_H!(it::Class_iterate, timer::class_advanced_timer, pars::Class_parameters)
    start_advanced_timer(timer, "CACHE/update_H")
    it.cache.H = eval_lag_hess(it.nlp, it.point.x, it.point.y, 1.0)
    it.cache.H_updated = true

    nzH = nonzeros(it.cache.H)
    if isbad(nzH)
        println("NaN or Inf")
        throw(Eval_NaN_error(getbad(nzH), it.point.x, "H"))
    end
    pause_advanced_timer(timer, "CACHE/update_H")
end

function update_grad!(it::Class_iterate, timer::class_advanced_timer, pars::Class_parameters)
    start_advanced_timer(timer, "CACHE/update_grad")
    it.cache.grad = eval_grad_f(it.nlp, it.point.x)
    it.cache.grad_updated = true
    if isbad(it.cache.grad)
        println("NaN or Inf")
        throw(Eval_NaN_error(getbad(it.cache.grad), it.point.x, "grad"))
    end
    pause_advanced_timer(timer, "CACHE/update_grad")
end

function update_J!(it::Class_iterate, timer::class_advanced_timer, pars::Class_parameters)
  start_advanced_timer(timer, "CACHE/update_J")
  @assert(length(it.point.x) == length(it.cache.grad))
  it.cache.J = eval_jac(it.nlp, it.point.x)
  it.cache.J_updated = true
  if isbad(nonzeros(it.cache.J))
      println("NaN or Inf")
      throw(Eval_NaN_error(getbad(it.cache.J), it.point.x, "J"))
  end
  pause_advanced_timer(timer, "CACHE/update_J")
end
