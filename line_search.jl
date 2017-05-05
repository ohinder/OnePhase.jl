type Class_ls_info
    step_size_P::Float64
    step_size_D::Float64
    num_steps::Int64

    function Class_ls_info()
        return new(NaN,NaN,0)
    end

    function Class_ls_info(step_size_P,step_size_D,num_steps)
        return new(step_size_P,step_size_D,num_steps)
    end
end

function scale_direction(dir::Class_point, step_size::Float64)
    new_dir = deepcopy(dir)
    shrink_direction!(new_dir, step_size)
    return new_dir
end

#=function max_step_LP(iter::Class_iterate, dir::Class_point, pars::Class_parameters)
    for i = 1:ncon(iter)

    end
end=#
function simple_max_step_LP(iter::Class_iterate, dir::Class_point, threshold::Float64)
    q = 1.0 / (1.0 - threshold)
    ratio = maximum( [1.0; -q * dir.s ./ iter.point.s; -q * dir.y ./ iter.point.y] );
    return 1.0 / ratio
end

function accept_aggressive(intial_it::Class_iterate, candidate::Class_iterate, dir::Class_point, pars::Class_parameters)
    if is_feasible(candidate, pars)
        return :success
    else
        return :infeasible
    end
end

function accept_stable(iter::Class_iterate, candidate::Class_iterate, dir::Class_point, pars::Class_parameters)
    if !is_feasible(candidate, pars)
        return :infeasible
    end

    predict_red = merit_function_predicted_reduction(iter, dir);
    if predict_red >= 0.0
        warn("predicted reduction non-negative")
        return :predict_red_non_negative
    end

    old_merit = eval_merit_function(iter)
    new_merit = eval_merit_function(candidate)
    actual_red = new_merit - old_merit

    #@show old_merit, new_merit

    if actual_red < predict_red * pars.predict_reduction_factor
        return :success
    else
        return :not_enough_progress
    end
end

function simple_ls(iter::Class_iterate, orginal_dir::Class_point, accept::Function, pars::Class_parameters)
    step_size = simple_max_step_LP(iter, orginal_dir, pars.fraction_to_boundary)

    for i = 1:pars.ls_num_backtracks
      dir = scale_direction(orginal_dir, step_size)
      candidate = move(iter, dir)
      status = accept(iter, candidate, dir, pars)

      if status != :not_enough_progress
          if status == :success
            return :success, candidate, Class_ls_info(step_size, step_size, i)
          else
            return status, iter, Class_ls_info(step_size, step_size, i)
          end
      end
      step_size *= pars.ls_backtracking_factor
    end

    return :max_ls_it, iter, Class_ls_info(step_size, step_size, pars.ls_num_backtracks)
end

function aggressive_ls(iter::Class_iterate, org_dir::Class_point, accept::Function, pars::Class_parameters)
    step_size = 0.5
    dir = scale_direction(org_dir, step_size)
    candidate = move(iter, dir)

    if accept(candidate, pars)
        step_size_list = 1.0 - [2.^(-linspace(2.0,10.0)); 0.0]
        for new_step_size in step_size_list
          dir = scale_direction(org_dir, new_step_size)
          new_candidate = move(iter, dir)

          if !accept(new_candidate, pars)
              return candidate, Class_ls_info(step_size,step_size)
          end
          step_size = new_step_size
        end

        return candidate, Class_ls_info(step_size,step_size)
    else
        for i = 1:30
          step_size *= pars.ls_backtracking_factor
          dir = scale_direction(org_dir, step_size)
          candidate = move(iter, dir)

          if accept(candidate, pars)
              return candidate, Class_ls_info(step_size,step_size)
          end
        end
    end

    @show step_size
    error("ls fails!")
end
