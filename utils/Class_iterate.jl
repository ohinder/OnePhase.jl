type Class_cache
    function Class_cache()
      return new()
    end
end

abstract abstract_nlp;

type Class_local_info
    delta::Float64
end

type Class_iterate
    point::Class_point
    nlp::abstract_nlp
    cache::Class_cache
    local_info::Class_local_info

    function Class_iterate(intial_point::Class_point, nlp::abstract_nlp, local_info::Class_local_info)
      return new(intial_point, nlp, Class_cache(), local_info)
    end
end

import Base.copy

function copy(it::Class_iterate)
   new_point = deepcopy(it.point)
   new_it = Class_iterate(new_point, it.nlp, it.local_info)
end

function validate(it::Class_iterate)
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


function move(it::Class_iterate, dir::Class_point, step_size::Float64, pars::Class_parameters)
    if pars.move_type == :primal_dual
      primal_res_old = (it.point.s - eval_a(it))
      new_it = copy(it)
      new_it.point.mu += dir.mu * step_size
      new_it.point.x += dir.x * step_size
      new_a = eval_a(new_it)
      new_it.point.s = new_a + (new_it.point.mu / it.point.mu) * primal_res_old
      #new_it.point.s += dir.s

      if minimum(new_it.point.s) > 0.0
        lb1, ub1 = dual_bounds(it, new_it.point.y, dir.y, pars.comp_feas)
        lb2, ub2 = dual_bounds(it, new_it.point.y, dir.y, pars.comp_feas_agg)

        if pars.move_primal_seperate_to_dual
            if lb1 < ub1 / 1.01
                ∇a = eval_jac(new_it)
                scale = dual_scale(new_it)

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
                  #step_size_D = (lb1 + ub1) / 2.0
                end

                step_size_D = sum(res .* q) / sum(q.^2)
                step_size_D = min(ub, max(lb, step_size_D))


                new_it.point.y += dir.y * step_size_D

                if !is_feasible(new_it, pars.comp_feas)
                  @show minimum(new_it.point.y), maximum(new_it.point.s)
                  println("line=",@__LINE__, "file=",@__FILE__)
                  warn("infeasibility should have been detected earlier!!!")
                  @show (lb1, ub1), (lb2, ub2)
                  return it, false, step_size_D
                else
                  return new_it, true, step_size_D
                end
            else
              return it, false, step_size
            end
        else
          if lb <= step_size && step_size <= ub
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
      new_it.point.s = new_a + (new_it.point.mu/it.point.mu) * (it.point.s - eval_a(it))
      #new_it.point.s += dir.s
      new_it.point.y = new_it.point.mu ./ new_it.point.s
    end
end


function move(it::Class_iterate, dir::Class_point, pars::Class_parameters)
    if pars.move_type == :primal_dual
      new_it = copy(it)
      new_it.point.mu += dir.mu
      new_it.point.x += dir.x
      new_a = eval_a(new_it)
      new_it.point.s = new_a + (new_it.point.mu/it.point.mu) * (it.point.s - eval_a(it))
      #new_it.point.s += dir.s
      new_it.point.y += dir.y
      #new_it.point.y = new_it.point.mu ./ new_it.point.s
    elseif pars.move_type == :primal
      new_it = copy(it)
      new_it.point.mu += dir.mu
      new_it.point.x += dir.x
      new_a = eval_a(new_it)
      new_it.point.s = new_a + (new_it.point.mu/it.point.mu) * (it.point.s - eval_a(it))
      #new_it.point.s += dir.s
      new_it.point.y = new_it.point.mu ./ new_it.point.s
    end

    return new_it;
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
