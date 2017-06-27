type Class_kkt_ls <: abstract_ls_info
    step_size_P::Float64
    step_size_D::Float64
    num_steps::Int64
    frac_progress::Float64
    predict_red::Float64
    actual_red::Float64
    cur_merit::Float64

    function Class_kkt_ls(iter::Class_iterate, dir::Class_point, pars::Class_parameters)
        this = new()
        this.predict_red = merit_function_predicted_reduction(iter, dir, 1.0);
        this.cur_merit = eval_merit_function(iter, pars)

        return this
    end

    function Class_kkt_ls(step_size_P::Float64,step_size_D::Float64,num_steps::Int64,frac_progress::Float64, predict_red::Float64)
        this = new(step_size_P, step_size_D, num_steps, frac_progress, predict_red)
        return this
    end
end

function accept_func_kkt!(accept::abstract_ls_info, iter::Class_iterate, candidate::Class_iterate, dir::Class_point, step_size::Float64, filter::Array{Class_filter,1}, pars::Class_parameters, timer::class_advanced_timer)

    if satisfies_filter!(filter, candidate, accept.step_size_P, pars, timer)
    #if kkt_err_cand / kkt_err < (1.0 - pars.kkt_reduction_factor * accept.step_size_P)
      return :success
    else
      return :not_enough_progress
    end
end


function accept_func!(accept::Class_kkt_ls, iter::Class_iterate, candidate::Class_iterate, dir::Class_point, step_size::Float64, filter::Array{Class_filter,1}, pars::Class_parameters, timer::class_advanced_timer)
    return accept_func_kkt!(accept, iter, candidate, dir, step_size, filter, pars, timer)
end
