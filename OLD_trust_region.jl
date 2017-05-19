function delta_with_correct_progress!(kkt_solver::abstract_KKT_system_solver, pars::Class_parameters)
    delta = kkt_solver.delta
    delta_lb = 0.0
    delta_ub = Inf

    for i = 1:10
        inertia = factor!(kkt_solver, delta)

        if inertia == 1
          compute_direction!(kkt_solver)

          if min_radius <= cur_info.radius && current_radius <= cur_info.radius
              return :success
          end
        else
          kkt_solver.lambda_lb = delta
          delta_lb = max(delta_lb, kkt_solver.lambda_lb)
        end
    end
end












type solve_attempt
    delta::Float64
    radius::Float64
    inertia::Bool
    progress_factor::Float64
end

function select_delta_via_regular_falsi(hist::Array{solve_attempt,1}, target_radius::Float64)
    if length(hist) > 1
      sort(hist, by = s -> s.radius) # sort by increasing radius
      @show hist

      if hist[1].radius > target_radius
        delta = hist[1].delta * 10.0
      elseif
        delta = hist[1].delta / 10.0
      else
        for i = 2:length(hist)
          if hist[i-1].radius <= target_radius && target_radius <= hist[i].radius
            if hist[i].inertia == 1
              delta = (hist[i].radius - hist[i-1].radius) / ( hist[i].delta - hist[i-1].delta );
            else
              # improve!
              delta = ( hist[i].delta + hist[i-1].delta ) / 2.0;
            end

            return delta
          end
        end
        error("regular falsi failed")
      end
    end

    error("regular falsi needs at least two points");
end

function trust_region_direction!(kkt_solver::abstract_KKT_system_solver, min_radius::Float64, max_radius::Float64, delta_lb::Float64, delta_ub::Float64, pars::Class_parameters)
    @assert(min_radius < max_radius)

    delta = kkt_solver.delta

    hist = Array{solve_attempt,1}()

    for i = 1:10
        target_radius = (min_radius + max_radius) / 2.0
        delta = select_delta(hist, target_radius)
        inertia = factor!(kkt_solver, delta)

        if inertia == 1
          compute_direction!(kkt_solver)
          cur_info = solve_attempt(delta, norm(kkt_solver.dir.x,2))
          push!(hist, cur_info);

          if min_radius <= cur_info.radius && current_radius <= cur_info.radius
              return :success
          end
        else
          kkt_solver.lambda_lb = delta
          delta_lb = max(delta_lb, kkt_solver.lambda_lb)
        end
    end
end
