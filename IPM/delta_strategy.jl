function ipopt_strategy!(iter::Class_iterate, kkt_solver::abstract_KKT_system_solver, pars::Class_parameters, timer::class_advanced_timer)
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
          println("delta=",rd(delta), "inertia=", pd(inertia))
        end

        if inertia == 1
          return :success, num_fac, delta
        elseif i == 1
          if get_delta(iter) != 0.0
            delta = max(DELTA_MIN, get_delta(iter) * pars.delta.dec)
          else
            delta = choose_delta_start(iter, kkt_solver, pars)
            #@show delta
            #delta = DELTA_START
          end
        else
            delta = delta * pars.delta.inc
        end

        if delta > DELTA_MAX
          dx = norm(kkt_solver.dir.x,Inf)
          dy = norm(kkt_solver.dir.y,Inf)
          ds = norm(kkt_solver.dir.s,Inf)

          my_warn("ipopt_strategy failed with delta_max=$DELTA_MAX, delta=$delta, i=$i")
          my_warn("num_fac=$num_fac, inertia=$inertia, status=$status, dir_x=$dx, dir_y=$dx, dir_s=$ds")

          #error("ipopt_strategy failed with too big a delta")
          return :failure, iter, delta
        end
    end

    error("max it")
end

function choose_delta_start(iter::Class_iterate, kkt_solver::abstract_KKT_system_solver, pars::Class_parameters)
    return pars.delta.start

    # do something clever to reduce number of factorizations

end
