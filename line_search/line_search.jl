abstract abstract_ls_info;

include("move.jl")
include("filter_ls.jl")
include("agg_ls.jl")
include("stable_ls.jl")
include("kkt_ls.jl")

function scale_direction(dir::Class_point, step_size::Float64)
    new_dir = deepcopy(dir)
    new_dir.x *= step_size
    new_dir.y *= step_size
    new_dir.s *= step_size
    new_dir.mu *= step_size
    new_dir.primal_scale *= step_size

    return new_dir
end


function simple_max_step(val::Array{Float64,1}, dir::Array{Float64,1}, lb::Array{Float64,1})
    ratio = maximum( [1.0; -dir ./ (val - lb) ] ) #; -q * dir.y ./ iter.point.y] );
    return 1.0 / ratio
end

#=
function max_step_primal_dual(iter::Class_iterate, dir::Class_point, threshold::Float64)
    return simple_max_step([iter.point.s; iter.point.y], [dir.s; dir.y], threshold)
end

function max_step_primal(iter::Class_iterate, dir::Class_point, threshold::Float64)
    return simple_max_step(iter.point.s, dir.s, threshold)
end

function simple_max_step(val::Array{Float64,1}, dir::Array{Float64,1}, threshold::Float64)
    q = 1.0 / (1.0 - threshold)
    ratio = maximum( [1.0; -q * dir ./ val ] ) #; -q * dir.y ./ iter.point.y] );
    return 1.0 / ratio
end=#

function Blank_ls_info()
    this = Class_stable_ls(0.0,0.0, 0, 0.0, 0.0)
    return this
end

#=function ls_feasible_solution(iter::Class_iterate, dir::Class_point, accept_type::Symbol, filter::Array{Class_filter,1},  pars::Class_parameters, min_step_size::Float64, timer::class_advanced_timer)

end=#

function lb_s_predict(iter,dir,pars)
  ex = pars.fraction_to_boundary_predict_exp
  x_thres = max(2.0 * norm(dir.x,Inf)^2, norm(dir.x,Inf)^ex)
  return min(pars.fraction_to_boundary_predict * iter.point.s, x_thres)
end

function basic_checks(iter::Class_iterate, f_it::Class_iterate, dir::Class_point)
  primal_accurate = norm(eval_jac_prod(f_it, dir.x) - eval_jac_prod(iter, dir.x),Inf) < norm(get_primal_res(iter),Inf) #* 10.0
  dual_accurate = norm(eval_jac_T_prod(f_it, dir.y) - eval_jac_T_prod(iter, dir.y),Inf) < norm(eval_grad_lag(iter,get_mu(iter)),Inf) #* 10.0
  #comp_accurate = -50.0 * get_mu(iter) < minimum(comp_predicted(iter,dir,1.0)) && maximum(comp_predicted(iter,dir,1.0)) < 50.0 * get_mu(iter)
  comp_accurate = norm(comp_predicted(iter,dir,1.0),Inf) < 50.0 * get_mu(iter)

  return comp_accurate #&& primal_accurate && dual_accurate
  #return true
end

function simple_ls(iter::Class_iterate, f_it::Class_iterate, dir::Class_point, accept_type::Symbol, filter::Array{Class_filter,1},  pars::Class_parameters, min_step_size::Float64, timer::class_advanced_timer)
    start_advanced_timer(timer, "SIMPLE_LS")

    # compute fraction to boundary
    lb_s = lb_s_predict(iter, dir, pars)
    step_size_P = simple_max_step(iter.point.s, dir.s, lb_s)

    if pars.max_step_primal_dual == true
      step_size_P = min(step_size_P,simple_max_step(iter.point.y, dir.y, 0.2 * iter.point.y))
      #step_size_P = max_step_primal(iter, dir, frac_to_bound)
    elseif pars.max_step_primal_dual != false
      error("SIMPLE_LS")
    end

    if accept_type == :accept_stable
      accept_obj = Class_stable_ls(iter, dir, pars)
    elseif accept_type == :accept_aggressive
      accept_obj = Class_agg_ls(iter, dir, pars)
    elseif accept_type == :accept_filter
      accept_obj = Class_filter_ls(iter, dir, pars)
    elseif accept_type == :accept_kkt
      accept_obj = Class_kkt_ls(iter, dir, pars)
    elseif accept_type == :accept_comp
      accept_obj = Class_comp_ls(iter, dir, pars)
    else
      error("acceptance function not defined")
    end

    accept_obj.num_steps = 0

    #if !basic_checks(iter, f_it, dir)
    #  return :basic_checks, iter, accept_obj
    #end

    if !accept_obj.do_ls
        if pars.LS_non_negative_predicted_gain
          return :predict_red_non_negative, iter, accept_obj
        else
          my_warn("predicted reduction non-negative")
        end
    end

    if pars.output_level >= 5
      println(pd("α_P"), pd("α_D"), pd("merit_diff"), pd("comp_diff"), pd("phi_diff"), pd("kkt_diff"), pd("dx"), pd("dy"), pd("ds"), pd("status"))
    end

    for i = 1:pars.ls_num_backtracks
      status = :none

      if step_size_P >= min_step_size
        start_advanced_timer(timer,"SIMPLE_LS/move/primal")
        candidate, move_status = move_primal(iter, dir, step_size_P, pars, timer)
        pause_advanced_timer(timer,"SIMPLE_LS/move/primal")

        accept_obj.step_size_P = step_size_P
        accept_obj.num_steps = i


        if move_status == :success
            start_advanced_timer(timer,"SIMPLE_LS/move/dual_bounds")
            lb, ub = dual_bounds(candidate, candidate.point.y, dir.y, pars.comp_feas)
            lb_y_new = min(dir.y .* (dir.s + norm(dir.x,Inf)^2) ./ candidate.point.s, pars.fraction_to_boundary * candidate.point.y)
            ub = min(ub, simple_max_step(candidate.point.y, dir.y, lb_y_new))
            pause_advanced_timer(timer,"SIMPLE_LS/move/dual_bounds")

            if lb >= ub
              move_status = :dual_infeasible
            end
            if !pars.move_primal_seperate_to_dual && !(lb <= step_size_P && step_size_P <= ub)
              move_status = :dual_infeasible
            end
        end

        if move_status == :success
            no_nan = update_J!(candidate, timer, pars) && update_grad!(candidate, timer, pars)
            if no_nan
              start_advanced_timer(timer,"SIMPLE_LS/move/dual")
              candidate, move_status, step_size_D = move_dual(candidate, dir, step_size_P, lb, ub, pars, timer)
              pause_advanced_timer(timer,"SIMPLE_LS/move/dual")
              accept_obj.step_size_D = step_size_D
            else
              move_status = :NaN_ERR
            end
        end


        if move_status == :success
          start_advanced_timer(timer,"SIMPLE_LS/accept?")
          no_nan = update_obj!(candidate, timer, pars)
          if no_nan
            status = accept_func!(accept_obj, iter, candidate, dir, step_size_P, filter, pars, timer)
          else
            status = :NaN_ERR
          end
          pause_advanced_timer(timer,"SIMPLE_LS/accept?")
        else
          status = move_status
        end

        @assert(is_updated_correction(iter.cache))

        # SOME DEBUGGING CODE
        if pars.debug_mode >= 1
          obj = copy(get_fval(iter))
          update_obj!(iter, timer, pars)
          #@show get_fval(iter), obj, get_fval(candidate)
          if(obj != get_fval(iter))
              error("this shouldn't happen!!!")
          end
        end

        # OUTPUT INFORMATION ABOUT LINE SEARCH
        if pars.output_level >= 5
          diff = eval_merit_function(candidate, pars) - eval_merit_function(iter, pars)
          comp_diff = norm(comp(candidate),Inf) - norm(comp(iter), Inf)
          mu = iter.point.mu
          phi_diff = eval_phi(candidate, mu) - eval_phi(iter, mu)

          dx = step_size_P * norm(dir.x,2);
          dy = norm(candidate.point.y - iter.point.y,2); ds = norm(candidate.point.s - iter.point.s,2);
          kkt_diff = norm(eval_grad_lag(candidate, get_mu(candidate)),Inf) / norm(eval_grad_lag(iter, get_mu(iter)),Inf)
          println(rd(step_size_P), rd(step_size_D), rd(diff), rd(comp_diff), rd(phi_diff), rd(kkt_diff), rd(dx), rd(dy), rd(ds), pd(status))
        end

        #if is_feas
        #  pause_advanced_timer(timer,"SIMPLE_LS/feas")
        #end

        if status == :success
          pause_advanced_timer(timer,"SIMPLE_LS")
          return :success, candidate, accept_obj
        elseif status == :predict_red_non_negative
          pause_advanced_timer(timer, "SIMPLE_LS")
          return status, iter, accept_obj
        end

        #if status == :s_bound
        #  step_size_P *= pars.ls_backtracking_factor # min(pars.ls_backtracking_factor * step_size_P, 1.0 / maximum( (get_s(iter) - get_s(candidate)) ./ get_s(iter) ))
        #else
        step_size_P *= pars.ls_backtracking_factor
        #end
      else
        if pars.output_level >= 5
          println(rd(step_size_P), pd("N/A"), pd("N/A"), pd(:min_α))
        end

        pause_advanced_timer(timer, "SIMPLE_LS")
        return :min_α, iter, accept_obj
      end
    end

    pause_advanced_timer(timer, "SIMPLE_LS")
    return :max_ls_it, iter, accept_obj
end



function eigenvector_ls(iter::Class_iterate, dir::Class_point, pars::Class_parameters)
    step_size_P = 1.0

    best_candidate = iter;
    intial_val = eval_merit_function(iter, pars)
    best_val = intial_val

    max_it = 10;

    i = 1;
    for i = 1:max_it
      candidate_pos, is_feas, step_size_D = move(iter, dir, step_size_P, pars, timer)
      candidate_neg, is_feas, step_size_D = move(iter, dir, -step_size_P, pars, timer)

      better = false

      new_val = eval_merit_function(candidate_pos, pars)
      if new_val < best_val
        best_candidate = candidate_pos
        best_val = new_val
        better = true
      end

      new_val = eval_merit_function(candidate_neg, pars)
      if new_val < best_val
        best_candidate = candidate_neg
        best_val = new_val
        better = true
      end

      if !better
        break
      end

      step_size_P *= 6.0
    end

    @show step_size_P, norm(dir.x), best_val - intial_val

    if i == max_it
      my_warn("max it reached for eig search")
    end

    return :success, best_candidate
end
