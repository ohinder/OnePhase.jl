abstract abstract_nlp;

type Class_cache
    fval::Float64
    cons::Array{Float64,1}
    grad::Array{Float64,1}
    J::SparseMatrixCSC{Float64,Int32}
    H::SparseMatrixCSC{Float64,Int32}

    function Class_cache()
      return new()
    end
end

type Class_local_info
    delta::Float64
    radius::Float64
end

type Class_iterate
    primal_residual_intial::Array{Float64,1}
    point::Class_point
    nlp::abstract_nlp
    cache::Class_cache
    local_info::Class_local_info
    ncon::Int64
    nvar::Int64

    function Class_iterate(intial_point::Class_point, nlp::abstract_nlp, local_info::Class_local_info, timer::class_advanced_timer)
      this = new(zeros(length(intial_point.x)), intial_point, nlp, Class_cache(), local_info)
      this.nvar = length(intial_point.x)
      this.ncon = length(intial_point.y)
      @assert(this.ncon == length(intial_point.y))

      this.cache = Class_cache()
      update_grad!(this, timer)
      update_cons!(this, timer)
      update_obj!(this, timer)
      update_H!(this, timer)
      update_J!(this, timer)

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

        ub_dyi = mu / (comp_feas * s) - y[i]
        lb_dyi = mu * comp_feas / s - y[i]

        if dy[i] > 0.0
          lb = max( lb_dyi / dy[i], lb)
          ub = min( ub_dyi / dy[i], ub)

          #if lb < ub
          #  @assert( y[i] + ub * dy[i] <= mu / (comp_feas * s) )
          #  @assert( y[i] + lb * dy[i] >= comp_feas * mu / s )
          #end
        elseif dy[i] < 0.0
          lb = max(ub_dyi / dy[i], lb)
          ub = min(lb_dyi / dy[i], ub)

          #if lb < ub
          #  @assert( y[i] + lb * dy[i] <= mu / (comp_feas * s) )
          #  @assert( y[i] + ub * dy[i] >= comp_feas * mu / s )
          #end
        elseif lb_dyi >= 0.0 || ub_dyi <= 0.0
          lb = 0.0
          ub = -1.0
        end
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
      update_cons!(new_it, timer)

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
                if lb2 < ub2 && false
                    scale = dual_scale(new_it, pars)
                    ∇a = get_jac(new_it)
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


                if lb2 < ub2 && false
                  step_size_D = ub2
                else
                  step_size_D = sum(res .* q) / sum(q.^2)
                  step_size_D = ub1 #min(ub1, max(lb1, step_size_D))
                end


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
    new_cache.J = cache.J
    new_cache.H = cache.H
    new_cache.cons = cache.cons
    new_cache.fval = cache.fval
    new_cache.grad = cache.grad

    it.cache = new_cache

    pause_advanced_timer(timer, "CACHE/copy")
end

function update_cons!(it::Class_iterate, timer::class_advanced_timer)
    start_advanced_timer(timer, "CACHE/update_cons")
    it.cache.cons = eval_a(it.nlp, it.point.x)
    pause_advanced_timer(timer, "CACHE/update_cons")
end

function update_obj!(it::Class_iterate, timer::class_advanced_timer)
    start_advanced_timer(timer, "CACHE/update_obj")
    it.cache.fval = eval_f(it.nlp, it.point.x)
    pause_advanced_timer(timer, "CACHE/update_obj")
end

function update_H!(it::Class_iterate, timer::class_advanced_timer)
    start_advanced_timer(timer, "CACHE/update_H")
    it.cache.H = eval_lag_hess(it.nlp, it.point.x, it.point.y, 1.0)
    pause_advanced_timer(timer, "CACHE/update_H")
end

function update_grad!(it::Class_iterate, timer::class_advanced_timer)
    start_advanced_timer(timer, "CACHE/update_grad")
    it.cache.grad = eval_grad_f(it.nlp, it.point.x)
    pause_advanced_timer(timer, "CACHE/update_grad")
end

function update_J!(it::Class_iterate, timer::class_advanced_timer)
  start_advanced_timer(timer, "CACHE/update_J")
  @assert(length(it.point.x) == length(it.cache.grad))
  it.cache.J = eval_jac(it.nlp, it.point.x)
  pause_advanced_timer(timer, "CACHE/update_J")
end
