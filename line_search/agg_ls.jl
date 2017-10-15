
type Class_agg_ls <: abstract_ls_info
    step_size_P::Float64
    step_size_D::Float64
    num_steps::Int64
    frac_progress::Float64
    predict_red::Float64
    actual_red::Float64

    function Class_agg_ls(iter::Class_iterate, dir::Class_point, pars::Class_parameters)
        this = new()

        this.step_size_D = 0.0
        this.step_size_P = 0.0
        this.num_steps = 0
        this.predict_red = -1.0
        this.frac_progress = NaN
        this.actual_red = NaN

        return this
    end
end


function accept_func!(accept::Class_agg_ls, intial_it::Class_iterate, candidate::Class_iterate, dir::Class_point, step_size::Float64, filter::Array{Class_filter,1},  pars::Class_parameters, timer::class_advanced_timer)
    if scaled_dual_feas(candidate, pars) < get_mu(candidate) * pars.agg_protect_factor
      return :success
    else
      return :not_enough_progress
    end
end
