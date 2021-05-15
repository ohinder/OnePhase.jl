
mutable struct Class_agg_ls <: abstract_ls_info
    step_size_P::Float64
    step_size_D::Float64
    num_steps::Int64
    frac_progress::Float64
    predict_red::Float64
    actual_red::Float64
    do_ls::Bool

    function Class_agg_ls(iter::Class_iterate, dir::Class_point, pars::Class_parameters)
        this = new()

        this.step_size_D = 0.0
        this.step_size_P = 0.0
        this.num_steps = 0
        this.predict_red = -1.0
        this.frac_progress = NaN
        this.actual_red = NaN
        #mu_ = get_mu(iter) - dir.mu
        #this.do_ls = 0.1 < minimum(comp_predicted_scaled(iter,dir,1.0)) && maximum(comp_predicted(iter,dir,1.0)) .< 10 * get_mu(iter)
        #comp_max = (get_mu(iter) - dir.mu) / pars.ls.comp_feas_agg + get_mu(iter) / 4.0
        #this.do_ls = true
        eta = - dir.mu / get_mu(iter)
        gamma = 1.0 - eta
        y_tilde = (gamma * get_mu(iter) - eta * iter.point.y .* get_primal_res(iter)) ./ iter.point.s
        this.do_ls = dot(eval_grad_lag(iter, get_mu(iter) * gamma, y_tilde), dir.x) < 0.0

        #this.do_ls = true
        #LinearAlgebra.norm(comp_predicted(iter,dir,1.0),Inf) < get_mu(iter) / pars.ls.comp_feas_agg
        return this
    end
end


function accept_func!(accept::Class_agg_ls, intial_it::Class_iterate, candidate::Class_iterate, dir::Class_point, step_size::Float64, filter::Array{Class_filter,1},  pars::Class_parameters, timer::class_advanced_timer)
    tau = get_mu(candidate) / (scaled_dual_feas(candidate, pars) * (1.0 - pars.agg_protection_factor))
    max_thres = pars.agg_protection_factor
    #@show get_mu(candidate), get_mu(intial_it)
    #@show get_mu(candidate), get_mu(intial_it)
    if get_mu(candidate) / get_mu(intial_it) >= 1.0 - max_thres || tau >= 1.0
      #@show get_mu(candidate) / get_mu(intial_it), tau
      return :success, step_size
    else
      #println("failure")
      suggested_step_size = max(max_thres^2, accept.step_size_P * tau^2)
      return :not_enough_progress, suggested_step_size
    end
end
