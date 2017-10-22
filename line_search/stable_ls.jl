
type Class_stable_ls <: abstract_ls_info
    step_size_P::Float64
    step_size_D::Float64
    num_steps::Int64
    predict_red::Float64
    frac_progress::Float64
    actual_red::Float64
    cur_merit::Float64
    do_ls::Bool

    function Class_stable_ls(iter::Class_iterate, dir::Class_point, pars::Class_parameters)
        this = new()

        this.step_size_D = 0.0
        this.step_size_P = 0.0
        this.num_steps = 0
        this.predict_red = merit_function_predicted_reduction(iter, dir, 1.0);
        this.frac_progress = NaN
        this.actual_red = NaN
        this.cur_merit = eval_merit_function(iter, pars)
        this.do_ls = this.predict_red >= 0.0

        return this
    end

    function Class_stable_ls(step_size_P::Float64,step_size_D::Float64,num_steps::Int64,frac_progress::Float64, predict_red::Float64)
        this = new(step_size_P, step_size_D, num_steps, frac_progress, predict_red)
        return this
    end
end


function accept_func_stable!(accept::abstract_ls_info, iter::Class_iterate, candidate::Class_iterate, dir::Class_point, step_size::Float64, pars::Class_parameters, timer::class_advanced_timer)
  if accept.predict_red < 0.0
    #old_merit = eval_merit_function(iter, pars)
    old_merit = accept.cur_merit
    new_merit = eval_merit_function(candidate, pars)

    accept.actual_red = new_merit - old_merit

    #@show accept.predict_red, accept.actual_red
    #predict_red = merit_function_predicted_reduction(iter, dir, step_size);
    predict_red = accept.predict_red * step_size

    if accept.actual_red > 0.0
      return :negative_progress
    else
      accept.frac_progress = accept.actual_red / predict_red
      if accept.frac_progress > pars.predict_reduction_factor
          #@show accept.frac_progress, accept.actual_red
          return :success
      else
          return :not_enough_progress
      end
    end
  else
      return :not_enough_progress
  end
end


function accept_func!(accept::Class_stable_ls, iter::Class_iterate, candidate::Class_iterate, dir::Class_point, step_size::Float64, filter::Array{Class_filter,1}, pars::Class_parameters, timer::class_advanced_timer)
    status = accept_func_stable!(accept, iter, candidate, dir, step_size, pars, timer)
    return status
end
