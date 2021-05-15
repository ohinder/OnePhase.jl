#using LinearAlgebra
#using SparseArrays
mutable struct Schur_KKT_solver <: abstract_schur_solver
    # abstract_KKT_system_solver
    ls_solver::abstract_linear_system_solver # the linear system solver we wish to use, see the folder `linear_system_solvers`
    factor_it::Class_iterate # iterate where the factorization is computed
    delta_x_vec::Array{Float64,1} # the amount we perturb the Hessian
    delta_s_vec::Array{Float64,1} # MAYBE DELETE? I DON'T REALLY USE THIS ANYMORE??????
    rhs::System_rhs # the right hand side
    dir::Class_point # the direction that was computed
    kkt_err_norm::Class_kkt_error
    rhs_norm::Float64
    pars::Class_parameters
    schur_diag::Array{Float64,1} # diagonal values of the schur complement, if any of these elements is negative matrix inertia is incorrect.

    ready::Symbol # is the KKT solver correctly intialized

    # Schur_KKT_solver only
    Q::SparseMatrixCSC{Float64,Int64}
    #M::SparseMatrixCSC{Float64,Int64} # Is two linear systems really necessary????
    #K::SparseMatrixCSC{Float64,Int64}
    current_it::Class_iterate # current iterate
    reduct_factors::Class_reduction_factors # the amount that we want to reduce

    function Schur_KKT_solver()
      this = new()
      this.ready = :not_ready

      return this
    end
end


function kkt_associate_rhs!(kkt_solver::abstract_schur_solver, iter::Class_iterate, reduct_factors::Class_reduction_factors, timer::class_advanced_timer)
    start_advanced_timer(timer, "KKT/rhs");

    kkt_solver.rhs = System_rhs(iter, reduct_factors)
    kkt_solver.dir.mu = -(1.0 - reduct_factors.mu) * get_mu(iter)
    kkt_solver.dir.primal_scale = -(1.0 - reduct_factors.P) * iter.point.primal_scale

    kkt_solver.reduct_factors = reduct_factors
    kkt_solver.current_it = iter

    pause_advanced_timer(timer, "KKT/rhs");
end

function form_system!(kkt_solver::abstract_schur_solver, iter::Class_iterate, timer::class_advanced_timer)
    # form the matrix before we factorize it

    start_advanced_timer(timer, "SCHUR")
    start_advanced_timer(timer, "SCHUR/form_system");

    # TODO build specialized schur complement code.
    kkt_solver.Q = eval_J_T_J(iter, iter.point.y ./ iter.point.s) + get_lag_hess(iter);
    kkt_solver.schur_diag = diag(kkt_solver.Q)
    kkt_solver.factor_it = iter;
    kkt_solver.ready = :system_formed

    pause_advanced_timer(timer, "SCHUR/form_system");
    pause_advanced_timer(timer, "SCHUR");
end

function update_delta_vecs!(kkt_solver::abstract_schur_solver, delta_x_vec::Array{Float64,1}, delta_s_vec::Array{Float64,1}, timer::class_advanced_timer)
    start_advanced_timer(timer, "SCHUR")
    start_advanced_timer(timer, "SCHUR/delta_vecs")
    kkt_solver.delta_x_vec = delta_x_vec
    kkt_solver.delta_s_vec = delta_s_vec

    if sum(abs.(delta_s_vec)) > 0.0
        error("Not implemented")
        #kkt_solver.Q[i,i] =
        #kkt_solver.M + eval_J_T_J(kkt_solver.factor_it, delta_s_vec) + spdiagm(delta_x_vec)
    else
        for i = 1:size(kkt_solver.Q,1)
            kkt_solver.Q[i,i] = kkt_solver.schur_diag[i] + delta_x_vec[i]
        end
    end

    kkt_solver.ready = :delta_updated
    pause_advanced_timer(timer, "SCHUR/delta_vecs")
    pause_advanced_timer(timer, "SCHUR")
end

function factor_implementation!(kkt_solver::abstract_schur_solver, timer::class_advanced_timer)
    return ls_factor!(kkt_solver.ls_solver, kkt_solver.Q, dim(kkt_solver.factor_it), 0, timer)
end

function compute_direction_implementation!(kkt_solver::Schur_KKT_solver, timer::class_advanced_timer)
    # compute a search direction
    start_advanced_timer(timer, "SCHUR")

    factor_it = kkt_solver.factor_it
    #∇a_org = get_jac(factor_it);
    y_org = get_y(factor_it);
    s_org = get_s(factor_it);

    rhs = kkt_solver.rhs

    #r1 + ∇a_org' * (( r3 + (r2 .* y_org) ) ./ s_org)
    start_advanced_timer(timer, "SCHUR/rhs");
    symmetric_primal_rhs = rhs.primal_r + rhs.comp_r ./ y_org
    Σ_vec = ( y_org ./ s_org )
    y_ = (rhs.primal_r .* Σ_vec + rhs.comp_r ./ s_org )
    schur_rhs = rhs.dual_r + eval_jac_T_prod(factor_it, y_)
    pause_advanced_timer(timer, "SCHUR/rhs");

    dir = kkt_solver.dir
    dir.x = solver_schur_rhs(schur_rhs, kkt_solver, timer)

    # there are two ways to update s and y
    if true
      #J = get_jac(kkt_solver.current_it)
      dir.y = -(eval_jac_prod(factor_it,dir.x) - symmetric_primal_rhs) .* Σ_vec
      #dir.y = -(J * dir.x - symmetric_primal_rhs) .* Σ_vec
      #dir.s = ( rhs.comp_r - dir.y .* s_org ) ./ y_org
      dir.s = eval_jac_prod(factor_it,dir.x) - rhs.primal_r
    else
      dir.s = eval_jac_prod(factor_it,dir.x) - rhs.primal_r
      dir.y = ( rhs.comp_r - dir.s .* y_org ) ./ s_org # (mu_target - s_cur .* y_cur)
    end

    check_for_nan(dir)

    start_advanced_timer(timer, "SCHUR/kkt_err");
    update_kkt_error!(kkt_solver, Inf, timer)
    pause_advanced_timer(timer, "SCHUR/kkt_err");
    pause_advanced_timer(timer, "SCHUR")
end


function solver_schur_rhs(schur_rhs::Vector, kkt_solver::abstract_schur_solver, timer::class_advanced_timer)
  # Solve Ax = b where b is formed by the function compute_direction_implementation!().
  # The main purpose of this function is to do iterative refinement.

  fit = kkt_solver.factor_it
  #∇a_org = get_jac(factor_it);
  y_org = get_y(fit);
  s_org = get_s(fit);
  #(symmetric_primal_rhs .* Σ_vec);Σ
  dir = kkt_solver.dir;
  pars = kkt_solver.pars

  Σ_vec = ( y_org ./ s_org )
  sqrt_Σ_vec = sqrt.( y_org ./ s_org )

  # generalize!!!
  output_level = pars.output_level
  res_old = schur_rhs
  if output_level >= 4
    println("res", 0, " ", rd(LinearAlgebra.norm(res_old,2)))
  end

  dir_x = zeros(length(dir.x));
  if pars.kkt.ItRefine_BigFloat
    dir_x = convert(Array{BigFloat,1}, dir_x)
  end

  for i = 1:pars.kkt.ItRefine_Num
      if output_level >= 4
        println("res", i, " ", rd(Float64(LinearAlgebra.norm(res_old,2))))
      end

      start_advanced_timer(timer, "SCHUR/iterative_refinement");
      dir_x += ls_solve(kkt_solver.ls_solver, res_old, timer);

      start_advanced_timer(timer, "SCHUR/iterative_refinement/residual");
      jac_res = eval_jac_T_prod( fit , Σ_vec .* eval_jac_prod(fit, dir_x) )
      hess_res = hess_product(fit, dir_x) + kkt_solver.delta_x_vec .* dir_x
      res_old = schur_rhs - ( jac_res + hess_res )
      pause_advanced_timer(timer, "SCHUR/iterative_refinement/residual");

      #res_old = res
      pause_advanced_timer(timer, "SCHUR/iterative_refinement");
  end

  if output_level >= 4
    println("res ", rd(Float64(LinearAlgebra.norm(res_old,2))))
  end

  return dir_x
end
