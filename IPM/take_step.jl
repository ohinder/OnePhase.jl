
function probe(iter::Class_iterate, kkt_solver::abstract_KKT_system_solver, pars::Class_parameters, timer::class_advanced_timer)
    kkt_associate_rhs!(kkt_solver, iter, Reduct_affine(), timer)
    compute_direction!(kkt_solver, timer)

    lb_s = lb_s_predict(iter, kkt_solver.dir, pars)

    m = length(iter.point.s)
    alpha_s = simple_max_step(iter.point.s, kkt_solver.dir.s, lb_s)
    alpha_y = simple_max_step(iter.point.y, kkt_solver.dir.y, zeros(m))
    return min(alpha_s, alpha_y)
end


function take_step!(iter::Class_iterate, eta::Class_reduction_factors, kkt_solver::abstract_KKT_system_solver, ls_mode::Symbol, filter::Array{Class_filter,1}, pars::Class_parameters, min_step_size::Float64, timer::class_advanced_timer)
    kkt_associate_rhs!(kkt_solver, iter, eta, timer)
    compute_direction!(kkt_solver, timer)

    if false #kkt_solver.kkt_err_norm.ratio > 1e-1
        @show kkt_solver.kkt_err_norm
        return :failure, iter, Class_ls_info()
    end

    success, iter, step_info = simple_ls(iter, kkt_solver.dir, ls_mode, filter, pars, min_step_size, timer)

    return success, iter, step_info
end

function take_step2!( be_aggressive::Bool, iter::Class_iterate, kkt_solver::abstract_KKT_system_solver, filter::Array{Class_filter,1}, pars::Class_parameters, timer::class_advanced_timer)
    if be_aggressive
        if pars.agg_eta == :mehrotra || pars.agg_eta == :mehrotra_stb
          σ = probe(iter, kkt_solver, pars, timer)
          #gamma = (1.0 - σ)
          gamma = min(0.5, (1.0 - σ)^2)
          #gamma = (1.0 - σ) #min(0.5, )
          #gamma = min(0.5,(1.0 - σ)^1.5)
          if pars.agg_eta == :mehrotra
            reduct_factors = Class_reduction_factors(gamma, gamma, gamma)
          elseif pars.agg_eta == :mehrotra_stb
            reduct_factors = Class_reduction_factors(gamma, 0.0, gamma)
          end
        else
          if pars.agg_eta == :affine
            reduct_factors = Reduct_affine()
          elseif pars.agg_eta == :constant
            reduct_factors = Class_reduction_factors(0.1, 0.1, 0.1)
          else
            error("pars.agg_eta choice not recognized")
          end
        end

        ls_mode = pars.ls_mode_agg
        min_step_size = pars.min_step_size_agg_ratio * min(1.0, 1.0 / maximum(- get_primal_res(iter) ./ iter.point.s))
    else
        ls_mode = pars.ls_mode_stable_correction
        reduct_factors = Reduct_stable()
        min_step_size = pars.min_step_size_stable
    end

    kkt_associate_rhs!(kkt_solver, iter, reduct_factors, timer)
    compute_direction!(kkt_solver, timer)

    if false #kkt_solver.kkt_err_norm.ratio > 1e-1
        @show kkt_solver.kkt_err_norm
        return :failure, iter, Class_ls_info(), reduct_factors
    end

    success, iter, step_info = simple_ls(iter, kkt_solver.factor_it, kkt_solver.dir, ls_mode, filter, pars, min_step_size, timer)
    return success, iter, step_info, reduct_factors
end
