type Schur_KKT_solver <: abstract_schur_solver
    # abstract_KKT_system_solver
    ls_solver::abstract_linear_system_solver
    factor_it::Class_iterate
    delta_x_vec::Array{Float64,1}
    delta_s_vec::Array{Float64,1}
    rhs::System_rhs
    dir::Class_point
    kkt_err_norm::Class_kkt_error
    rhs_norm::Float64
    pars::Class_parameters

    ready::Symbol #

    # Schur_KKT_solver only
    true_diag::Array{Float64,1}
    M::SparseMatrixCSC{Float64,Int64}
    K::SparseMatrixCSC{Float64,Int64}
    current_it::Class_iterate
    reduct_factors::Class_reduction_factors

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
    start_advanced_timer(timer, "SCHUR")
    start_advanced_timer(timer, "SCHUR/form_system");
    x = iter.point.x;
    y = iter.point.y
    s = iter.point.s

    ## REMEMBER  M is triangular!!!
    # this could be sped up significantly!!!
    kkt_solver.M = eval_J_T_J(iter, iter.point.y ./ iter.point.s) + get_lag_hess(iter);
    kkt_solver.true_diag = diag(kkt_solver.M)
    kkt_solver.factor_it = iter;
    kkt_solver.ready = :system_formed

    pause_advanced_timer(timer, "SCHUR/form_system");
    pause_advanced_timer(timer, "SCHUR");
end

function update_diag!(Mat::SparseMatrixCSC{Float64,Int64})

end

function update_delta_vecs!(kkt_solver::abstract_schur_solver, delta_x_vec::Array{Float64,1}, delta_s_vec::Array{Float64,1}, timer::class_advanced_timer)
    start_advanced_timer(timer, "SCHUR")
    start_advanced_timer(timer, "SCHUR/delta_vecs")
    kkt_solver.delta_x_vec = delta_x_vec
    kkt_solver.delta_s_vec = delta_s_vec

    if sum(abs(delta_s_vec)) > 0.0
        kkt_solver.K = kkt_solver.M + eval_J_T_J(kkt_solver.factor_it, delta_s_vec) + spdiagm(delta_x_vec)
    else
        kkt_solver.K = kkt_solver.M + spdiagm(delta_x_vec)
    end

    kkt_solver.ready = :delta_updated
    pause_advanced_timer(timer, "SCHUR/delta_vecs")
    pause_advanced_timer(timer, "SCHUR")
end

function factor_implementation!(kkt_solver::abstract_schur_solver, timer::class_advanced_timer)
    return ls_factor!(kkt_solver.ls_solver, kkt_solver.K, dim(kkt_solver.factor_it), 0, timer)
end

function compute_direction_implementation!(kkt_solver::Schur_KKT_solver, timer::class_advanced_timer)
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
  fit = kkt_solver.factor_it
  #∇a_org = get_jac(factor_it);
  y_org = get_y(fit);
  s_org = get_s(fit);
  #(symmetric_primal_rhs .* Σ_vec);Σ
  dir = kkt_solver.dir;
  pars = kkt_solver.pars

  Σ_vec = ( y_org ./ s_org )

  # generalize!!!
  output_level = pars.output_level
  res_old = schur_rhs
  if output_level >= 4
    println("res", 0, " ", rd(norm(res_old,2)))
  end

  dir_x = zeros(length(dir.x));
  if pars.ItRefine_BigFloat
    dir_x = convert(Array{BigFloat,1}, dir_x)
  end

  for i = 1:pars.ItRefine_Num
      start_advanced_timer(timer, "SCHUR/iterative_refinement");
      dir_x += ls_solve(kkt_solver.ls_solver, res_old, timer)[:];

      start_advanced_timer(timer, "SCHUR/iterative_refinement/residual");
      jac_res = eval_jac_T_prod( fit , Σ_vec .* eval_jac_prod(fit, dir_x) )
      hess_res = hess_product(fit, dir_x) + kkt_solver.delta_x_vec .* dir_x
      res = schur_rhs - ( jac_res + hess_res )
      pause_advanced_timer(timer, "SCHUR/iterative_refinement/residual");

      if output_level >= 4
        println("res", i, " ", rd(Float64(norm(res,2))))
      end

      res_old = res
      pause_advanced_timer(timer, "SCHUR/iterative_refinement");
  end

  return dir_x
end
