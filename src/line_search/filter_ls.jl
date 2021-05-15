
#function merit_function_anticipated_reduction(iter, dir)
#    return norm(comp(iter),Inf)^3 / (iter.point.mu)^2 + norm(dir.x,2)^2 * get_delta(iter)
#end

mutable struct Class_filter_ls <: abstract_ls_info
    step_size_P::Float64
    step_size_D::Float64
    num_steps::Int64
    frac_progress::Float64
    predict_red::Float64
    actual_red::Float64
    cur_merit::Float64
    do_ls::Bool

    function Class_filter_ls(iter::Class_iterate, dir::Class_point, pars::Class_parameters)
        this = new()
        this.step_size_D = 0.0
        this.step_size_P = 0.0
        this.num_steps = 0

        g = eval_grad_phi(iter, get_mu(iter))
        this.predict_red = -comp_merit(iter) + 0.5 * (dot(g,dir.x) - get_delta(iter) * norm(dir.x,2)^2)
        #merit_function_predicted_reduction(iter, dir, 1.0);
        this.cur_merit = eval_merit_function(iter, pars)

        sufficient_descent = dot(g,dir.x) < 0.0 #-0.5 * norm(g,2)^2 / norm(get_lag_hess(iter),2)^2
        merit_reduce = merit_function_predicted_reduction(iter, dir, 1.0) < this.predict_red / 2.0
        #comp_sufficient = -get_mu(iter) < minimum(comp_predicted(iter,dir,1.0)) && maximum(comp_predicted(iter,dir,1.0)) < 100.0 * get_mu(iter)
        #comp_sufficient = -get_mu(iter) * pars.ls.comp_feas_agg < minimum(comp_predicted(iter,dir,1.0)) && maximum(comp_predicted(iter,dir,1.0)) < get_mu(iter) * (1.0 / pars.ls.comp_feas_agg - 1.0)
        this.do_ls = sufficient_descent #&& merit_reduce
        #&& comp_sufficient
        #-0.5 * merit_function_anticipated_reduction(iter, dir)

        return this
    end

    function Class_filter_ls(step_size_P::Float64,step_size_D::Float64,num_steps::Int64,frac_progress::Float64, predict_red::Float64)
        this = new(step_size_P, step_size_D, num_steps, frac_progress, predict_red)
        return this
    end
end

mutable struct Class_filter
    fval::Float64
    scaled_kkt_err::Float64
    mu::Float64

    function Class_filter(iter::Class_iterate, pars::Class_parameters)
      this = new()
      this.fval = eval_merit_function(iter, pars)

      kkt_err = norm(eval_grad_lag(iter, get_mu(iter)),Inf)
      if pars.ls.kkt_include_comp
         kkt_err += norm(comp(iter),Inf)
      end
      this.scaled_kkt_err = kkt_err * dual_scale(iter, pars) #eval_kkt_err(iter, pars)
      this.mu = iter.point.primal_scale

      return this
    end
end

function add!(ar::Array{Class_filter,1}, iter::Class_iterate, pars::Class_parameters )
    push!( ar, Class_filter(iter, pars) )
end

function satisfies_filter!(ar::Array{Class_filter,1}, can::Class_iterate, step_size::Float64, pars::Class_parameters, timer::class_advanced_timer)
    p = Class_filter(can, pars)

    failure = false

    for i = 1:length(ar)
        #tol = 1e-3 * max(1.0, p.scaled_kkt_err)
        expected_reduction = pars.ls.kkt_reduction_factor
        #
        kkt_reduction = (p.scaled_kkt_err / ar[i].scaled_kkt_err < (1.0 - pars.ls.kkt_reduction_factor * step_size))
        fval_reduction = p.fval < ar[i].fval - (p.scaled_kkt_err)^2
        fval_no_increase = p.fval < ar[i].fval + sqrt(p.scaled_kkt_err) #+ p.scaled_kkt_err^2
        net_reduction = p.scaled_kkt_err + p.fval < ar[i].fval + ar[i].scaled_kkt_err - (p.scaled_kkt_err)^2

        if pars.ls.filter_type == :default
            failure = !(p.mu < ar[i].mu || kkt_reduction)
        elseif pars.ls.filter_type == :test1
            failure = !(p.mu < ar[i].mu || kkt_reduction || fval_reduction)
        elseif pars.ls.filter_type == :test2
            failure = !(p.mu < ar[i].mu || (kkt_reduction && fval_no_increase)) # || fval_reduction)
        elseif pars.ls.filter_type == :test3
            failure = !(p.mu < ar[i].mu || net_reduction)
        else
            error("unknown filter type!!!")
        end

        if failure #&& p.fval < ar[i].fval + (p.scaled_kkt_err) ))
          # && p.fval < ar[i].fval + (p.scaled_kkt_err)))
            #println(p.scaled_kkt_err / ar[i].scaled_kkt_err)
            # p.fval < ar[i].fval - (p.scaled_kkt_err)^2 ||
            return false
        end
    end

    return true
end



function accept_func!(accept::Class_filter_ls, iter::Class_iterate, candidate::Class_iterate, dir::Class_point, step_size::Float64, filter::Array{Class_filter,1}, pars::Class_parameters, timer::class_advanced_timer)
    status, suggested_step_size_P = accept_func_stable!(accept, iter, candidate, dir, step_size, pars, timer)

    if status == :success || status == :infeasible
      return status, suggested_step_size_P
    else
      return accept_func_kkt!(accept, iter, candidate, dir, step_size, filter, pars, timer)
    end
end
