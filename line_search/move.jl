

function move_primal(it::Class_iterate, dir::Class_point, step_size::Float64, pars::Class_parameters, timer::class_advanced_timer)
    if pars.move_type == :primal_dual
      new_it = copy(it, timer)
      new_it.point.x += dir.x * step_size

      if !update_cons!(new_it, timer, pars)
        return new_it, :NaN_ERR
      end

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

      lb_s = min(it.point.s * pars.fraction_to_boundary, norm(dir.x,Inf)^2)
      if !all(new_it.point.s .>= lb_s)
          return new_it, :s_bound
      end

      new_it.point.mu += dir.mu * step_size

      return new_it, :success
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
    end

    #boundary = min(0.01 * y, abs(dy) .* abs(dy) ./ it.point.s
    #ub = min(ub, simple_max_step(y, dy, 10.0^(-3.0) * y))

    return lb, ub
end

function move_dual(new_it::Class_iterate, dir::Class_point, step_size_P::Float64, lb::Float64, ub::Float64, pars::Class_parameters, timer::class_advanced_timer)
    if pars.move_primal_seperate_to_dual
      small_step = max(lb,min(ub,step_size_P))
      if pars.dual_ls == 2
        scale = dual_scale(new_it, pars)
        old_y = deepcopy(new_it.point.y)
        intial_value = eval_kkt_err(new_it, pars)
        new_it.point.y = old_y + small_step * dir.y
        small_value = eval_kkt_err(new_it, pars)
        new_it.point.y = old_y + ub * dir.y
        big_value = eval_kkt_err(new_it, pars)

        if big_value < intial_value * (1.0 - pars.kkt_reduction_factor)
          step_size_D = ub
        else
          step_size_D = small_step
        end

        new_it.point.y = old_y
      elseif pars.dual_ls == 1 || pars.dual_ls == 3
        scale_D = dual_scale(new_it, pars)
        scale_mu = dual_scale(new_it, pars)
        #scale = 1.0;

        #∇a = get_jac(new_it)
        dual_res = eval_grad_lag(new_it, get_mu(new_it))
        q = [scale_D * eval_jac_T_prod(new_it,dir.y); scale_mu * new_it.point.s .* dir.y];
        prox_dual_res = dual_res + get_delta(new_it) * dir.x * step_size_P * (pars.dual_ls == 3)
        #@show get_delta(new_it)
        res = [scale_D * dual_res; -scale_mu * comp(new_it)]

        step_size_D = sum(res .* q) / sum(q.^2)

        step_size_D = min(step_size_D,ub)
        step_size_D = max(step_size_D,small_step)
      else
        step_size_D = ub
      end
    else
      step_size_D = step_size_P
    end

    new_it.point.y += dir.y * step_size_D

    if !is_feasible(new_it, pars.comp_feas)
      @show minimum(new_it.point.y), maximum(new_it.point.s)
      println("line=",@__LINE__, "file=",@__FILE__)
      my_warn("infeasibility should have been detected earlier!!!")
      @show step_size_D
      @show (lb, ub)
      return new_it, :dual_infeasible, step_size_D
    else
      return new_it, :success, step_size_D
    end
end
