type Class_reduction_factors
    P::Float64
    D::Float64
    mu::Float64
    function Class_reduction_factors(P::Float64,D::Float64,mu::Float64)
      return new(P,D,mu)
    end
    function Class_reduction_factors()
      return new(NaN,NaN,NaN)
    end
end

type System_rhs
    dual_r::Array{Float64,1}
    primal_r::Array{Float64,1}
    comp_r::Array{Float64,1}
    function System_rhs()
        return new()
    end

    function System_rhs(it::Class_iterate, reduct::Class_reduction_factors)
      dual_target = -eval_grad_lag(it) * (1.0 - reduct.D)
      primal_target = -eval_primal_residual(it) * (1.0 - reduct.P)
      mu_target = get_mu(it) * reduct.mu
      s = get_s(it)
      y = get_y(it)

      return new(dual_target, primal_target, mu_target - s .* y);
    end
end

import Base.LinAlg.norm
function norm(rhs::System_rhs, p::Float64)
    return norm([rhs.dual_r, rhs.primal_r, rhs.comp_r], p)
end

abstract abstract_KKT_system_solver;

function initialize!(kkt_solver::abstract_KKT_system_solver, intial_it::Class_iterate)
    initialize!(kkt_solver.ls_solver)
    kkt_solver.dir = zero_point(dim(intial_it),ncon(intial_it))
end

function predicted_lag_change(fact_it::Class_iterate, delta::Float64, dir::Class_point)
    return delta * dir.x + eval_lag_hess(fact_it) * dir.x - eval_jac(fact_it)' * dir.y
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

function update_kkt_error!(ss::abstract_KKT_system_solver, iter::Class_iterate, p::Float64)
    factor_iter = ss.factor_it;
    rhs = ss.rhs;
    dir = ss.dir;

    candidate = move(iter,ss.dir)
    dual_scaling_can = dual_scale(candidate)

    error_D = (predicted_lag_change(factor_iter, ss.delta, dir) - rhs.dual_r) * dual_scaling_can;
    error_P = eval_jac(factor_iter) * dir.x - dir.s - rhs.primal_r;
    error_mu = get_s(factor_iter) .* dir.y + get_y(factor_iter) .* dir.s  - rhs.comp_r;

    #@show norm(error_D), norm(error_P), norm(error_mu)
    overall = norm([error_D; error_P; error_mu], p)

    dual_scaling_org = dual_scale(iter)
    rhs_norm = norm([rhs.dual_r * dual_scaling_org; rhs.primal_r; rhs.comp_r], p)

    ratio = overall / rhs_norm;

    ss.kkt_err_norm = Class_kkt_error(norm(error_D,p), norm(error_P,p), norm(error_mu,p), overall, rhs_norm, ratio)
end

function factor_delta!()
    # move factor! to here.
end

function factor!(kkt_solver::abstract_KKT_system_solver, iter::Class_iterate)
    start_advanced_timer("kkt/factor");
    # is this the best strategy ? what if delta is dynamic? ? ? ?
    form_system!(kkt_solver, iter)

    max_it = 100

    delta = get_delta(iter)
    i = 0
    for i = 1:max_it
      kkt_solver.delta = delta
      inertia = factor!(kkt_solver, delta)

      if inertia == 1
        if i == 1
          set_delta(iter, delta / 3.0)
        end

        break
      else
        #if i == 1
        #  set_delta(iter, max(1e-8, get_delta(iter) / 3.0))
        #else
        set_delta(iter, max(1e-8, get_delta(iter) * 10.0) )
        #end
        delta = get_delta(iter)
      end
    end

    if i == max_it
      @show delta
      error("delta too large!")
    end
    pause_advanced_timer("kkt/factor");
end

#function compute_direction!(kkt_solver::abstract_KKT_system_solver, rhs::System_rhs, dir_mu::Float64)
#
#end

function compute_direction!(kkt_solver::abstract_KKT_system_solver, iter::Class_iterate, eta::Class_reduction_factors)
    start_advanced_timer("kkt/direction");

    mu_dir = -(1.0 - eta.mu) * get_mu(iter)
    kkt_solver.rhs = System_rhs(iter, eta)
    kkt_solver.dir.mu = mu_dir

    compute_direction!(kkt_solver)
    #kkt_solver.rhs_norm = norm(rhs, Inf)

    update_kkt_error!(kkt_solver, iter, Inf)

    pause_advanced_timer("kkt/direction");
end

function pick_KKT_solver(pars::Class_parameters)
  kkt_solver_type = pars.kkt_solver_type
  linear_solver_type = pars.linear_solver_type
  safe = pars.linear_solver_safe_mode
  if kkt_solver_type == :symmetric
    my_kkt_solver = Symmetric_KKT_solver()
    if linear_solver_type == :julia
      my_kkt_solver.ls_solver = linear_solver_JULIA(:unsymmetric, safe)
    elseif linear_solver_type == :mumps
      my_kkt_solver.ls_solver = linear_solver_MUMPS(:symmetric, safe)
    end
  elseif kkt_solver_type == :schur
    my_kkt_solver = Schur_KKT_solver()
    if linear_solver_type == :julia
      my_kkt_solver.ls_solver = linear_solver_JULIA(:definite, safe)
    elseif linear_solver_type == :mumps
      my_kkt_solver.ls_solver = linear_solver_MUMPS(:definite, safe)
    end
  else
      error("pick a solver!")
  end

  my_kkt_solver.kkt_err_norm = Class_kkt_error()

  return my_kkt_solver
end


include("schur.jl")
include("symmetric.jl")

#function compute_direction(kkt_solver::KKT_system_solver, iter::Class_iterate, eta::Float64)
#
#end
