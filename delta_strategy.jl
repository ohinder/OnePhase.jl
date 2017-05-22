function ipopt_strategy!(iter::Class_iterate, stable_reduct_factors::Class_reduction_factors, kkt_solver::abstract_KKT_system_solver, ls_mode::Symbol, filter::Array{Class_filter,1}, pars::Class_parameters)

    form_system!(kkt_solver, iter)
    kkt_associate_rhs!(kkt_solver, iter, stable_reduct_factors)

    step_info = Blank_ls_info()
    new_iter = nothing

    MAX_IT = 30
    DELTA_MIN = 1e-14
    DELTA_START = 1e-6

    delta = 0.0

    if pars.output_level >= 4
        println(pd("delta"), pd("inertia"), pd("step_size_P"))
    end

    for i = 1:MAX_IT
        success = false

        if pars.use_delta_s
          inertia = factor!( kkt_solver, delta , delta )
        else
          inertia = factor!( kkt_solver, delta )
        end

        if inertia == 1
          compute_direction!(kkt_solver)

          if kkt_solver.kkt_err_norm.ratio < pars.saddle_err_tol
            #if delta == 0.0
            #  status, new_iter, step_info = simple_ls(iter, kkt_solver.dir, :accept_kkt, filter, pars, pars.min_step_size_stable)
            #else
            status, new_iter, step_info = simple_ls(iter, kkt_solver.dir, ls_mode, filter, pars, pars.min_step_size_stable)
            #end

            if status == :success
              success = true
            end
          end
        end

        if pars.output_level >= 4
          println(rd(delta), pd(inertia), rd(step_info.step_size_P))
        end

        if success && step_info.step_size_P > pars.min_step_size_stable
          #((step_info.step_size_P > 1e-2 && KKT_THRESHOLD <= delta) || (step_info.step_size_P > 1e-4 && delta < KKT_THRESHOLD ))
            set_delta(new_iter, delta)
            return status, new_iter, step_info
        elseif i == 1
          if get_delta(iter) != 0.0
            delta = max(DELTA_MIN, get_delta(iter) / pi)
          else
            delta = DELTA_START
          end
        else
            delta = delta * 8.0
        end
    end

    my_warn("ipopt_strategy failed with delta=$delta")
    return :failure, iter, Blank_ls_info()
end



function optimal_stable_step!(iter::Class_iterate, stable_reduct_factors::Class_reduction_factors, kkt_solver::abstract_KKT_system_solver, pars::Class_parameters)

    factor_at_approx_min_eigenvalue!(kkt_solver, iter)
    associate_rhs!(kkt_solver, iter, stable_reduct_factors)

    i = 1;

    max_it = 20
    step_info = nothing
    new_iter = nothing

    for i = 1:max_it
        if i == 2
          println("large delta regime")
        end
        compute_direction!(kkt_solver)
        status, new_iter, step_info = simple_ls(iter, kkt_solver.dir, :accept_stable, pars)

        if status == :success && step_info.step_size_P > 0.3
            break
        else
            set_delta( iter, get_delta(iter) * 2.0 )
            factor!( kkt_solver, get_delta(iter) )
        end
    end

    if i > max_it
        error("delta too big")
    end

    mode = :trust_region

    if i == 1 && get_delta(iter) > 1e-4 && mode == :trust_region
        println("trust region")
        alpha = 0.8

        old_candidate = new_iter
        step_info_old_candidate = step_info
        i = 1;
        for i = 1:max_it
          lambda = compute_eigenvector!(kkt_solver, iter)

          delta_candidate = (get_delta(iter) * (1 - alpha) - lambda * alpha)
          inertia = factor!( kkt_solver, delta_candidate )

          if inertia == 1
            compute_direction!(kkt_solver)
            @show norm(kkt_solver.dir.x)
            status, candidate, step_info_candidate = simple_ls(iter, kkt_solver.dir, :accept_stable, pars)
          else
            my_warn("eigenvector inaccurate")
          end

          if inertia == 0 || step_info_candidate.step_size_P < 1.0
              iter = old_candidate
              step_info = step_info_old_candidate
              break
          end
          old_candidate = candidate
          step_info_old_candidate = step_info_candidate
          set_delta( iter, delta_candidate )
        end

        if i == max_it
          my_warn("delta too small")
        end
    elseif get_delta(iter) > 1e-4 && mode == :eigenvector
        println("eigenvector")
        iter = new_iter
        dir_x_norm = norm(kkt_solver.dir.x,2)
        lambda = compute_eigenvector!(kkt_solver, iter)
        @show lambda
        shrink_direction!(kkt_solver.dir, norm(dir_x_norm,2) / 10.0)
        status, iter = eigenvector_ls(iter, kkt_solver.dir, pars)

        @show eigmin(full(Symmetric(kkt_solver.M)))

        H = eval_phi_hess(iter)
        @show eigmin(full(Symmetric(H)))
    else
        iter = new_iter
    end




    return :success, iter, step_info

end



function optimal_stable_step_v2!(iter::Class_iterate, stable_reduct_factors::Class_reduction_factors, kkt_solver::abstract_KKT_system_solver, pars::Class_parameters)
    form_system!(kkt_solver, iter)
    associate_rhs!(kkt_solver, iter, stable_reduct_factors)

    inertia = factor!( kkt_solver, 0.0 )
    if inertia == 1
        status, new_iter, step_info = simple_ls(iter, kkt_solver.dir, :accept_stable, pars)
        if status == :success
            set_delta(new_iter, 0.0 )
            return status, new_iter, step_info
        end
    end

    MAX_IT = 20
    DELTA_MIN = 1e-8
    DELTA_START = 100.0

    if get_delta(iter) == 0.0
      delta = DELTA_START
    else
      delta = max(get_delta(iter), DELTA_MIN)
    end

    lambda_lb = 0.0
    delta_lb = 0.0
    delta_ub = Inf

    increase_delta = false

    println("")
    println(pd("status"), pd("frac"), pd("prog"), pd("|dx|"), pd("delta"), pd("delta_lb"), pd("delta_ub"), pd("lambda_lb"))

    for i = 1:MAX_IT
        inertia = factor!( kkt_solver, delta )
        if inertia == 1
            compute_direction!(kkt_solver)
            status, new_iter, step_info = simple_ls(iter, kkt_solver.dir, :accept_stable, pars)

            #@show
            println(pd(status), rd(step_info.frac_progress), rd(step_info.actual_red),  rd(norm(kkt_solver.dir.x,2)), rd(delta), rd(delta_lb), rd(delta_ub), rd(lambda_lb))

            if status == :success
              if delta - lambda_lb > 1e-8 &&   pars.predict_reduction_factor_MAX < step_info.frac_progress #&& step_info.frac_progress < 1.0/pars.predict_reduction_factor_MAX

                delta_ub = delta
                increase_delta = false
              elseif step_info.frac_progress < pars.predict_reduction_factor
                delta_lb = delta
                increase_delta = true
              else
                set_delta(new_iter, delta)
                return status, new_iter, step_info # ACCEPT THE STEP
              end
            else
              delta_lb = delta
              increase_delta = true
            end
        else
            lambda_lb = delta
            delta_lb = delta
            increase_delta = true
        end

        @assert(lambda_lb <= delta_lb)
        @assert(delta_lb < delta_ub)

        if (delta_ub / delta_lb) < 1.01
          set_delta(new_iter, delta)
          return status, new_iter, step_info # ACCEPT THE STEP
        end

        if increase_delta
          delta *= 10.0
        else
          delta /= 9.0
        end

        if delta <= delta_lb || delta >= delta_ub
            delta = (delta_ub + delta_lb) / 2.0
        end

        #@show delta_lb, delta_ub, delta
    end

    error("max it for delta algorithm reached")
end
