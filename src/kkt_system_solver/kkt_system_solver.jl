@compat abstract type abstract_KKT_system_solver end
@compat abstract type abstract_schur_solver <: abstract_KKT_system_solver end



function initialize!(kkt_solver::abstract_KKT_system_solver, intial_it::Class_iterate)
    initialize!(kkt_solver.ls_solver)
    kkt_solver.dir = zero_point(dim(intial_it),ncon(intial_it))
end

function predicted_lag_change(kkt_solver::abstract_KKT_system_solver)
    #J = get_jac(kkt_solver.factor_it)
    fi = kkt_solver.factor_it

    dir_x = kkt_solver.dir.x
    dir_y = kkt_solver.dir.y

    tmp1 = eval_jac_prod(fi, dir_x)
    tmp2 = kkt_solver.delta_s_vec .* tmp1
    tmp3 = eval_jac_T_prod(fi, tmp2)
    delta_err = kkt_solver.delta_x_vec .* dir_x + tmp3
    J_err = eval_jac_T_prod(fi, dir_y)
    H_err = hess_product(kkt_solver.factor_it, dir_x)
    return delta_err + H_err - J_err
end

type Class_kkt_error
    error_D::Float64
    error_P::Float64
    error_mu::Float64
    overall::Float64
    rhs_norm::Float64
    ratio::Float64

    function Class_kkt_error(error_D::Float64,error_P::Float64,error_mu::Float64,overall::Float64, rhs_norm::Float64, ratio::Float64)
          return new(error_D,error_P,error_mu,overall,rhs_norm, ratio)
    end

    function Class_kkt_error()
        return new(0.0,0.0,0.0,0.0)
    end
end

function update_kkt_error!(ss::abstract_KKT_system_solver, p::Float64, timer::class_advanced_timer)
    factor_iter = ss.factor_it;
    rhs = ss.rhs;
    dir = ss.dir;

    dual_scaling_can = 1.0 #dual_scale(candidate)

    start_advanced_timer(timer, "update_kkt_error/compute/D")
    error_D = (predicted_lag_change(ss) - rhs.dual_r) * dual_scaling_can;
    pause_advanced_timer(timer, "update_kkt_error/compute/D")

    start_advanced_timer(timer, "update_kkt_error/compute/P")
    error_P = eval_jac_prod(factor_iter,dir.x) - dir.s - rhs.primal_r;
    pause_advanced_timer(timer, "update_kkt_error/compute/P")

    start_advanced_timer(timer, "update_kkt_error/compute/mu")
    error_mu = get_s(factor_iter) .* dir.y + get_y(factor_iter) .* dir.s  - rhs.comp_r;
    pause_advanced_timer(timer, "update_kkt_error/compute/mu")

    #@show norm(error_D), norm(error_P), norm(error_mu)
    overall = norm([error_D; error_P; error_mu], p)

    dual_scaling_org = 1.0 #dual_scale(iter)
    rhs_norm = norm([rhs.dual_r * dual_scaling_org; rhs.primal_r; rhs.comp_r], p)

    ratio = overall / rhs_norm;
    #@show ratio

    ss.kkt_err_norm = Class_kkt_error(norm(error_D,p), norm(error_P,p), norm(error_mu,p), overall, rhs_norm, ratio)
end

function factor!(kkt_solver::abstract_KKT_system_solver, delta_x::Float64, delta_s::Float64, timer::class_advanced_timer)
    update_delta!(kkt_solver, delta_x, delta_s, timer)
    factor!(kkt_solver, timer)
end

function factor!(kkt_solver::abstract_KKT_system_solver, delta_x::Float64, timer::class_advanced_timer)
    factor!(kkt_solver, delta_x, 0.0, timer)
end

function update_delta!(kkt_solver::abstract_KKT_system_solver, delta_x::Float64, delta_s::Float64, timer::class_advanced_timer)
  delta_x_vec = delta_x * ones(dim(kkt_solver.factor_it))
  delta_s_vec = delta_s * get_s(kkt_solver.factor_it).^(-2)
  update_delta_vecs!(kkt_solver, delta_x_vec, delta_s_vec, timer)
end

function factor_at_approx_min_eigenvalue!(kkt_solver::abstract_KKT_system_solver, iter::Class_iterate)
    start_advanced_timer("KKT/factor");
    form_system!(kkt_solver, iter)

    max_it = 100

    delta_min = 1e-8

    inertia = factor!(kkt_solver, 0.0)


    j = 1;
    if inertia == 1
      set_delta(iter, 0.0)
    else
      set_delta(iter, max(delta_min, get_delta(iter)))
      for i = 1:max_it
        inertia = factor!(kkt_solver, get_delta(iter))

        if inertia == 1
            set_delta(iter, get_delta(iter) / 20.0 )
        else
          break
        end
      end

      for j = 1:max_it
        inertia = factor!(kkt_solver, get_delta(iter))

        if inertia == 1
            break
        else
           set_delta(iter, get_delta(iter) * 3.0 )
        end
      end
    end

    if j == max_it
      @show get_delta(iter)
      error("delta too large!")
    end

    pause_advanced_timer("KKT/factor");
end


function kkt_associate_rhs!(kkt_solver::abstract_KKT_system_solver, iter::Class_iterate, eta::Class_reduction_factors, timer::class_advanced_timer)
    start_advanced_timer(timer, "KKT/rhs");

    kkt_solver.rhs = System_rhs(iter, eta)
    kkt_solver.dir.mu = -(1.0 - eta.mu) * get_mu(iter)
    kkt_solver.dir.primal_scale = -(1.0 - eta.P) * iter.point.primal_scale

    pause_advanced_timer(timer, "KKT/rhs");
end

function compute_direction!(kkt_solver::abstract_KKT_system_solver, timer::class_advanced_timer)
    start_advanced_timer(timer, "KKT/compute_direction");
    if kkt_solver.ready != :factored
        error("kkt solver not ready to compute direction!")
    end

    compute_direction_implementation!(kkt_solver, timer)
    check_for_nan(kkt_solver.dir)
    pause_advanced_timer(timer, "KKT/compute_direction");
end

function factor!(kkt_solver::abstract_KKT_system_solver, timer::class_advanced_timer)
    start_advanced_timer(timer, "KKT/factor");

    if kkt_solver.ready != :delta_updated
        error("kkt solver not ready to factor!")
    else
        kkt_solver.ready = :factored
    end
    inertia =  factor_implementation!(kkt_solver, timer)
    pause_advanced_timer(timer, "KKT/factor");

    return inertia
end

function compute_eigenvector!(kkt_solver::abstract_KKT_system_solver, iter::Class_iterate, timer::class_advanced_timer)
    start_advanced_timer(timer, "KKT/direction/eig");

    approx_eigvec = randn(dim(iter));

    for i = 1:20
      rhs = System_rhs(iter)
      rhs.dual_r = approx_eigvec
      kkt_solver.rhs = rhs
      kkt_solver.dir.mu = 0.0

      compute_direction!(kkt_solver, timer)
      approx_eigvec = kkt_solver.dir.x / norm(kkt_solver.dir.x,2)
    end

    lambda = norm(kkt_solver.dir.x,2)
    eig_est = (1 - kkt_solver.delta_x_vec[1] * lambda) / lambda

    kkt_solver.dir = scale_direction(kkt_solver.dir, 1.0 / lambda)

    pause_advanced_timer(timer, "KKT/direction/eig");

    return eig_est
end

function pick_KKT_solver(pars::Class_parameters)
  kkt_solver_type = pars.kkt.kkt_solver_type
  linear_solver_type = pars.kkt.linear_solver_type
  safe = pars.kkt.linear_solver_safe_mode
  recycle = pars.kkt.linear_solver_recycle

  if kkt_solver_type == :symmetric
    my_kkt_solver = Symmetric_KKT_solver()
    if linear_solver_type == :julia
      my_kkt_solver.ls_solver = linear_solver_JULIA(:unsymmetric, safe, recycle)
    elseif linear_solver_type == :mumps
      my_kkt_solver.ls_solver = linear_solver_MUMPS(:symmetric, safe, recycle)
    end
  elseif kkt_solver_type == :schur
    my_kkt_solver = Schur_KKT_solver()
    if linear_solver_type == :julia
      my_kkt_solver.ls_solver = linear_solver_JULIA(:definite, safe, recycle)
    elseif linear_solver_type == :mumps
      my_kkt_solver.ls_solver = linear_solver_MUMPS(:definite, safe, recycle)
    end
  elseif kkt_solver_type == :schur_direct
    my_kkt_solver = Schur_KKT_solver_direct()
    if linear_solver_type == :julia
      my_kkt_solver.ls_solver = linear_solver_JULIA(:definite, safe, recycle)
    elseif linear_solver_type == :mumps
      my_kkt_solver.ls_solver = linear_solver_MUMPS(:definite, safe, recycle)
    end
  else
      error("pick a solver!")
  end

  my_kkt_solver.kkt_err_norm = Class_kkt_error()
  my_kkt_solver.pars = pars

  return my_kkt_solver
end
