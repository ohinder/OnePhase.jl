function ipopt_strategy!(iter::Class_iterate, kkt_solver::abstract_KKT_system_solver, pars::Class_parameters, timer::class_advanced_timer)
    form_system!(kkt_solver, iter, timer)

    step_info = Blank_ls_info()

    MAX_IT = pars.delta.max_it
    DELTA_ZERO = pars.delta.zero
    DELTA_MIN = pars.delta.min
    DELTA_MAX = pars.delta.max

    delta = DELTA_ZERO

    if pars.output_level >= 4
        println(pd("delta"), pd("inertia"))
    end
    inertia = 0
    num_fac = 0
    status = :none

    for i = 1:MAX_IT
        if pars.use_delta_s
          inertia = factor!( kkt_solver, delta , delta, timer )
        else
          inertia = factor!( kkt_solver, delta, timer )
        end
        num_fac += 1

        if pars.output_level >= 4
          println(rd(delta), pd(inertia))
        end

        if inertia == 1
          set_delta(iter, delta)
          return :success, num_fac
        elseif i == 1
          if get_delta(iter) != 0.0
            delta = max(DELTA_MIN, get_delta(iter) / pi)
          else
            delta = choose_delta_start(iter, kkt_solver, pars)
            #@show delta
            #delta = DELTA_START
          end
        else
            delta = delta * 8.0
        end

        if delta > DELTA_MAX
          dx = norm(kkt_solver.dir.x,Inf)
          dy = norm(kkt_solver.dir.y,Inf)
          ds = norm(kkt_solver.dir.s,Inf)

          my_warn("ipopt_strategy failed with delta_max=$DELTA_MAX, delta=$delta, i=$i")
          my_warn("num_fac=$num_fac, inertia=$inertia, status=$status, dir_x=$dx, dir_y=$dx, dir_s=$ds")

          error("ipopt_strategy failed with too big a delta")
          return :failure, iter, Blank_ls_info()
        end
    end

    error("max it")
end

function choose_delta_start(iter::Class_iterate, kkt_solver::abstract_KKT_system_solver, pars::Class_parameters)
    return pars.delta.start

    # compute

end


#=
function ipopt_strategy_stable!(iter::Class_iterate, stable_reduct_factors::Class_reduction_factors, kkt_solver::abstract_KKT_system_solver, filter::Array{Class_filter,1}, pars::Class_parameters, timer::class_advanced_timer)
    form_system!(kkt_solver, iter, timer)
    kkt_associate_rhs!(kkt_solver, iter, stable_reduct_factors, timer)

    step_info = Blank_ls_info()
    new_iter = nothing

    MAX_IT = pars.delta.max_it
    DELTA_START = pars.delta.start
    DELTA_ZERO = pars.delta.zero #get_mu(iter) / ((1e2 + norm(get_x(iter),Inf)) * 1e2)
    DELTA_MIN = pars.delta.min
    DELTA_MAX = pars.delta.max #DELTA_MAX = 10.0 * norm(eval_lag_hess(iter), 1)

    delta = DELTA_ZERO

    if pars.output_level >= 4
        println(pd("delta"), pd("inertia"), pd("step_size_P"))
    end
    inertia = 0
    status = :none

    for i = 1:MAX_IT
        success = false

        if pars.use_delta_s
          inertia = factor!( kkt_solver, delta , delta, timer )
        else
          inertia = factor!( kkt_solver, delta, timer )
        end

        if inertia == 1
          compute_direction!(kkt_solver, timer)

          if kkt_solver.kkt_err_norm.ratio < pars.saddle_err_tol
            if delta <= DELTA_ZERO
              ls_mode = pars.ls_mode_stable_delta_zero
            else
              ls_mode = pars.ls_mode_stable_trust
            end

            status, new_iter, step_info = simple_ls(iter, kkt_solver.dir, ls_mode, filter, pars, pars.min_step_size_stable, timer)

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

        if DELTA_MAX < delta
          break
        end
    end

    dx = norm(kkt_solver.dir.x,Inf)
    dy = norm(kkt_solver.dir.y,Inf)
    ds = norm(kkt_solver.dir.s,Inf)

    my_warn("ipopt_strategy failed with delta_max=$DELTA_MAX, delta=$delta, inertia=$inertia, status=$status, dir_x=$dx, dir_y=$dx, dir_s=$ds")

    error("ipopt_strategy failed with to big a delta")
    return :failure, iter, Blank_ls_info()
end

function solve_with_radius(iter::Class_iterate, target_radius::Float64, kkt_solver::abstract_KKT_system_solver, pars::Class_parameters, timer::class_advanced_timer)
    new_iter = nothing

    MAX_IT = 50
    LAMBDA_LB = 0.0
    LAMBDA_UB = Inf

    DELTA = 1.0
    DELTA_LB = 0.0
    DELTA_UB = Inf

    radius = NaN
    best_dir = nothing
    println(pd("DELTA"), pd("DELTA_LB"), pd("DELTA_UB"), pd("RADIUS"))

    for i = 1:MAX_IT
        success = false

        if pars.use_delta_s
          inertia = factor!( kkt_solver, DELTA, DELTA, timer )
        else
          inertia = factor!( kkt_solver, DELTA, timer )
        end

        if inertia == 1
          compute_direction!(kkt_solver, timer)
          dir = kkt_solver.dir
          radius = norm(dir.x,2)
          if radius < target_radius / 2.0
            DELTA_UB = DELTA
          elseif radius > target_radius * 2.0
            DELTA_LB = DELTA
          else
            success = true
          end

          best_dir = dir
        else
          DELTA_LB = DELTA
          LAMBDA_LB = DELTA
          radius = NaN
        end

        println(rd(DELTA), rd(DELTA_LB), rd(DELTA_UB), rd(radius))

        set_delta(iter, DELTA)

        if success
          return true, best_dir
        end

        if DELTA_UB - DELTA_LB < 1e-8
          println("close delta")
          return false, best_dir
        end

        TRY_DELTA1 = (20.0 * DELTA_LB + 1e-8)
        TRY_DELTA2 = (DELTA_UB + DELTA_LB) / 2.0

        if TRY_DELTA1 <= TRY_DELTA2
          DELTA = TRY_DELTA1
        else
          DELTA = TRY_DELTA2
        end
    end

    println("trust region step failed!")
    return best_dir
end


function trust_region_strategy!(iter::Class_iterate, kkt_solver::abstract_KKT_system_solver, filter::Array{Class_filter,1}, pars::Class_parameters, timer::class_advanced_timer)
    form_system!(kkt_solver, iter, timer)
    kkt_associate_rhs!(kkt_solver, iter, Reduct_stable(), timer)


    DELTA_ZERO = 0.0
    inertia = factor!( kkt_solver, DELTA_ZERO, timer )

    if inertia == 1
        compute_direction!(kkt_solver, timer)
        status, new_iter, step_info = simple_ls(iter, kkt_solver.dir, pars.ls_mode_stable_delta_zero, filter, pars, pars.min_step_size_stable, timer)

        if status == :success
          set_radius( new_iter, norm(kkt_solver.dir.x,2) * step_info.step_size_P )
          set_prev_step( new_iter, :delta_zero )
          set_delta(new_iter, DELTA_ZERO)

          return status, new_iter, step_info
        end
    end

    if get_prev_step( iter ) == :init
        set_radius(iter, 1.0)
    end

    trust_status, dir = solve_with_radius(iter, get_radius(iter), kkt_solver, pars, timer)

    status, new_iter, step_info = simple_ls(iter, dir, pars.ls_mode_stable_trust, filter, pars, pars.min_step_size_stable, timer)

    if status == :success
        @show step_info.step_size_P, step_info.step_size_D
        if false
          if 0.1 < step_info.gain && step_info.gain < 0.95

          elseif step_info.gain >= 0.95
              set_radius(new_iter, get_radius(iter) * 3.0)
          elseif step_info.gain < 0.1
              set_radius(new_iter, get_radius(iter) / 4.0)
          end
        end

        if true
          if step_info.step_size_P == 1.0 && trust_status == true
            set_radius(new_iter, get_radius(iter) * 3.0)
          elseif step_info.step_size_P < 1.0
            set_radius(new_iter, get_radius(iter) * sqrt(step_info.step_size_P))
          end
        end

        set_prev_step( new_iter, :trust )

        @show get_radius(new_iter)
        @show get_prev_step(new_iter)

        return status, new_iter, step_info
    else
        set_radius(iter, get_radius(iter) * pars.min_step_size_stable )
        return status, iter, step_info
    end
end
=#
