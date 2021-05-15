#using LinearAlgebra
################################################################################
## Defines KKT system solvers
## KKT system solvers modify the orginal linear system to make it possible to use a linear system solver.
## For example,
## - schur.jl does a primal schur complement which can then be solved with a cholseky linear system solver.
## - symmetric.jl forms a symmetric system which can be solved using LDL
################################################################################

@compat abstract type abstract_KKT_system_solver
    # list variables
end
# The following functions should be defined for any KKT_solver <: abstract_KKT_system_solver:
# - form_system!(kkt_solver::abstract_schur_solver, iter::Class_iterate, timer::class_advanced_timer)
# - factor_implementation!(kkt_solver::KKT_solver, timer::class_advanced_timer)
# - compute_direction_implementation!(kkt_solver::KKT_solver, timer::class_advanced_timer)
# - update_delta_vecs!(kkt_solver::KKT_solver, delta_x_vec::Array{Float64,1}, delta_s_vec::Array{Float64,1}, timer::class_advanced_timer)

@compat abstract type abstract_schur_solver <: abstract_KKT_system_solver end

function initialize!(kkt_solver::abstract_KKT_system_solver, intial_it::Class_iterate)
    # call this before running using the kkt_solver
    initialize!(kkt_solver.ls_solver)
    kkt_solver.dir = zero_point(dim(intial_it),ncon(intial_it))
end

function predicted_lag_change(kkt_solver::abstract_KKT_system_solver)
    #################################################################################
    # INPUT:
    # kkt_solver with a direction and delta_x_vec representing how we perturbed the system to make the inertia correct.
    # OUTPUT:
    # vector denoting how do we expect gradient of the Lagragian to change
    #################################################################################

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

mutable struct Class_kkt_error
    # when we compute the directions how much error is there in the linear system.
    error_D::Float64 # norm of error in dual feasibility
    error_P::Float64 # norm error of error in primal feasibility
    error_mu::Float64 # norm error of error in complementarity
    overall::Float64 # total error
    rhs_norm::Float64 # norm of rhs of linear system
    ratio::Float64 # ratio of total error to norm of rhs. If this is less than zero then the direction is improving things.

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

    #@show LinearAlgebra.norm(error_D), LinearAlgebra.norm(error_P), LinearAlgebra.norm(error_mu)
    overall = LinearAlgebra.norm([error_D; error_P; error_mu], p)

    dual_scaling_org = 1.0 #dual_scale(iter)
    rhs_norm = LinearAlgebra.norm([rhs.dual_r * dual_scaling_org; rhs.primal_r; rhs.comp_r], p)

    ratio = overall / rhs_norm;
    #@show ratio

    ss.kkt_err_norm = Class_kkt_error(LinearAlgebra.norm(error_D,p), LinearAlgebra.norm(error_P,p), LinearAlgebra.norm(error_mu,p), overall, rhs_norm, ratio)
end

function factor!(kkt_solver::abstract_KKT_system_solver, delta_x::Float64, delta_s::Float64, timer::class_advanced_timer)
    # I don't think this is used anymore
    update_delta!(kkt_solver, delta_x, delta_s, timer)
    factor!(kkt_solver, timer)
end

function factor!(kkt_solver::abstract_KKT_system_solver, delta_x::Float64, timer::class_advanced_timer)
    # factorize the linear system for a particular perturbation to the Hessian of I * delta_x.
    factor!(kkt_solver, delta_x, 0.0, timer)
end

function update_delta!(kkt_solver::abstract_KKT_system_solver, delta_x::Float64, delta_s::Float64, timer::class_advanced_timer)
  delta_x_vec = delta_x * ones(dim(kkt_solver.factor_it))
  delta_s_vec = delta_s * get_s(kkt_solver.factor_it).^(-2)
  update_delta_vecs!(kkt_solver, delta_x_vec, delta_s_vec, timer)
end

function factor_at_approx_min_eigenvalue!(kkt_solver::abstract_KKT_system_solver, iter::Class_iterate)
    # Delete: I don't think this is used

    start_advanced_timer("KKT/factor");
    form_system!(kkt_solver, iter)

    max_it = 100

    delta_min = 1e-8

    inertia = factor!(kkt_solver, 0.0)


    #j = 1;
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
      
      counter_j = 1
      for j = 1:max_it
	counter_j = j
        
        inertia = factor!(kkt_solver, get_delta(iter))
        
        if inertia == 1
            break
        else
           set_delta(iter, get_delta(iter) * 3.0 )
        end

      end
    end

    if counter_j == max_it
      @show get_delta(iter)
      error("delta too large!")
    end

    pause_advanced_timer("KKT/factor");
end


function kkt_associate_rhs!(kkt_solver::abstract_KKT_system_solver, iter::Class_iterate, eta::Class_reduction_factors, timer::class_advanced_timer)
    # create right hand side for linear system and attach it to iter.
    start_advanced_timer(timer, "KKT/rhs");

    kkt_solver.rhs = System_rhs(iter, eta)
    kkt_solver.dir.mu = -(1.0 - eta.mu) * get_mu(iter)
    kkt_solver.dir.primal_scale = -(1.0 - eta.P) * iter.point.primal_scale

    pause_advanced_timer(timer, "KKT/rhs");
end

function compute_direction!(kkt_solver::abstract_KKT_system_solver, timer::class_advanced_timer)
    # compute a direction and store it in the kkt_solver object
    start_advanced_timer(timer, "KKT/compute_direction");
    if kkt_solver.ready != :factored
        error("kkt solver not ready to compute direction!")
    end

    compute_direction_implementation!(kkt_solver, timer)
    check_for_nan(kkt_solver.dir)
    pause_advanced_timer(timer, "KKT/compute_direction");
end

function factor!(kkt_solver::abstract_KKT_system_solver, timer::class_advanced_timer)
    # factorize the linear system
    # stores the factorize inside the kkt_solver
    start_advanced_timer(timer, "KKT/factor");

    if kkt_solver.ready != :delta_updated
        error("kkt solver not ready to factor kkt_solver.ready = $(kkt_solver.ready) != :delta_updated")
    else
        kkt_solver.ready = :factored
    end
    inertia = factor_implementation!(kkt_solver, timer)
    pause_advanced_timer(timer, "KKT/factor");

    return inertia # is the inertia correct?
end

function compute_eigenvector!(kkt_solver::abstract_KKT_system_solver, iter::Class_iterate, timer::class_advanced_timer)
    # Delete: I think this is old code not really used.
    start_advanced_timer(timer, "KKT/direction/eig");

    approx_eigvec = randn(dim(iter));

    for i = 1:20
      rhs = System_rhs(iter)
      rhs.dual_r = approx_eigvec
      kkt_solver.rhs = rhs
      kkt_solver.dir.mu = 0.0

      compute_direction!(kkt_solver, timer)
      approx_eigvec = kkt_solver.dir.x / LinearAlgebra.norm(kkt_solver.dir.x,2)
    end

    lambda = LinearAlgebra.norm(kkt_solver.dir.x,2)
    eig_est = (1 - kkt_solver.delta_x_vec[1] * lambda) / lambda

    kkt_solver.dir = scale_direction(kkt_solver.dir, 1.0 / lambda)

    pause_advanced_timer(timer, "KKT/direction/eig");

    return eig_est
end

function pick_KKT_solver(pars::Class_parameters)
    # selects the KKT_solver based on the parameters
  kkt_solver_type = pars.kkt.kkt_solver_type
  linear_solver_type = pars.kkt.linear_solver_type
  safe = pars.kkt.linear_solver_safe_mode
  recycle = pars.kkt.linear_solver_recycle

  if kkt_solver_type == :symmetric
        my_kkt_solver = Symmetric_KKT_solver()
        if linear_solver_type == :julia
          my_kkt_solver.ls_solver = linear_solver_JULIA(:symmetric, safe, recycle)
        elseif linear_solver_type == :mumps
          my_kkt_solver.ls_solver = linear_solver_MUMPS(:symmetric, safe, recycle)
        elseif linear_solver_type == :HSL
            my_kkt_solver.ls_solver = linear_solver_HSL(:symmetric, safe, recycle, pars.kkt.ma97_u)
        else
            error("pick a valid solver!")
        end
  elseif kkt_solver_type == :clever_symmetric
      my_kkt_solver = Clever_Symmetric_KKT_solver()
      if linear_solver_type == :julia
        my_kkt_solver.ls_solver = linear_solver_JULIA(:symmetric, safe, recycle)
      elseif linear_solver_type == :mumps
        my_kkt_solver.ls_solver = linear_solver_MUMPS(:symmetric, safe, recycle)
      elseif linear_solver_type == :HSL
          my_kkt_solver.ls_solver = linear_solver_HSL(:symmetric, safe, recycle, pars.kkt.ma97_u)
      else
          error("pick a valid solver!")
      end
  elseif kkt_solver_type == :schur
    my_kkt_solver = Schur_KKT_solver()
    if linear_solver_type == :julia
      my_kkt_solver.ls_solver = linear_solver_JULIA(:definite, safe, recycle)
    elseif linear_solver_type == :mumps
      my_kkt_solver.ls_solver = linear_solver_MUMPS(:definite, safe, recycle)
    else
      error("pick a valid solver!")
    end
  elseif kkt_solver_type == :schur_direct
    my_kkt_solver = Schur_KKT_solver_direct()
    if linear_solver_type == :julia
      my_kkt_solver.ls_solver = linear_solver_JULIA(:definite, safe, recycle)
    elseif linear_solver_type == :mumps
      my_kkt_solver.ls_solver = linear_solver_MUMPS(:definite, safe, recycle)
  else
    error("pick a valid solver!")
  end
  else
      error("pick a solver!")
  end

  my_kkt_solver.kkt_err_norm = Class_kkt_error()
  my_kkt_solver.pars = pars

  return my_kkt_solver
end

## helper functions

function diag_min(kkt_solver::abstract_KKT_system_solver)
    # find the smallest element of kkt_solver.schur_diag
    return minimum(kkt_solver.schur_diag)
end

function compute_schur_diag(iter::Class_iterate)
    # compute the diagonal of the primal schur complement matrix
    schur_scaling = get_y(iter) ./ get_s(iter);
    return diag(get_lag_hess(iter)) + eval_diag_J_T_J(iter,schur_scaling)
end
