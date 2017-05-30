type Schur_KKT_solver <: abstract_KKT_system_solver
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

    # Schur_KKT_solver only
    true_diag::Array{Float64,1}
    M::SparseMatrixCSC{Float64,Int64}
    K::SparseMatrixCSC{Float64,Int64}

    function Schur_KKT_solver()
      return new()
    end
end


function form_system!(kkt_solver::Schur_KKT_solver, iter::Class_iterate)
    start_advanced_timer("SCHUR/form_system");
    x = iter.point.x;
    y = iter.point.y
    s = iter.point.s

    ∇a = eval_jac(iter)

    #scaling = zeros(dim(iter))
    #for i = 1:dim(iter)
    #    scaling[i] = norm(∇a[:,i], Inf)
    #end

    #@show size(∇a)
    #@show maximum(scaling), minimum(scaling)
    #tic()

    Σ = spdiagm(iter.point.y ./ iter.point.s)

    kkt_solver.M = (∇a' * Σ * ∇a) + eval_lag_hess(iter);
    kkt_solver.true_diag = diag(kkt_solver.M)
    kkt_solver.factor_it = iter;
    pause_advanced_timer("SCHUR/form_system");

end



function update_delta_vecs!(kkt_solver::Schur_KKT_solver, delta_x_vec::Array{Float64,1}, delta_s_vec::Array{Float64,1})
    kkt_solver.delta_x_vec = delta_x_vec
    kkt_solver.delta_s_vec = delta_s_vec
    ∇a = ∇a = eval_jac(kkt_solver.factor_it)

    if sum(abs(delta_s_vec)) > 0.0
        kkt_solver.K = kkt_solver.M + (∇a' * spdiagm(delta_s_vec) * ∇a) + spdiagm(delta_x_vec)
    else
        kkt_solver.K = kkt_solver.M + spdiagm(delta_x_vec)
    end
end

function factor!(kkt_solver::Schur_KKT_solver)
    return ls_factor!(kkt_solver.ls_solver, kkt_solver.K, dim(kkt_solver.factor_it), 0)
end

function compute_direction!(kkt_solver::Schur_KKT_solver)
    factor_it = kkt_solver.factor_it
    ∇a_org = eval_jac(factor_it);
    H_org = eval_lag_hess(factor_it)
    y_org = get_y(factor_it);
    s_org = get_s(factor_it);

    rhs = kkt_solver.rhs

    #r1 + ∇a_org' * (( r3 + (r2 .* y_org) ) ./ s_org)
    symmetric_primal_rhs = rhs.primal_r + rhs.comp_r ./ y_org
    Σ_vec = ( y_org ./ s_org )
    schur_rhs = rhs.dual_r + ∇a_org' * (rhs.primal_r .* Σ_vec + rhs.comp_r ./ s_org )

    #(symmetric_primal_rhs .* Σ_vec);Σ
    dir = kkt_solver.dir;

    # generalize!!!
    output_level = kkt_solver.pars.output_level
    res_old = schur_rhs
    if output_level >= 4
      println("res", 0, " ", rd(norm(res_old,2)))
    end
    dir_x = convert(Array{BigFloat,1}, zeros(length(dir.x)))

    for i = 1:kkt_solver.pars.num_iterative_refinements
        dir_x += ls_solve(kkt_solver.ls_solver, res_old)[:];

        jac_res = ∇a_org' * ( Σ_vec  .* (∇a_org * dir_x) )
        res = schur_rhs - ( jac_res + H_org * dir_x + kkt_solver.delta_x_vec .* dir_x )

        if output_level >= 4
          println("res", i, " ", rd(Float64(norm(res,2))))
        end

        #if norm(res,2) > 0.95 * norm(res_old,2)
        #  break
        #end
        res_old = res
    end
    dir.x = dir_x

    if true
      dir.y = -(∇a_org * dir.x - symmetric_primal_rhs) .* Σ_vec
      dir.s = ( rhs.comp_r - dir.y .* s_org ) ./ y_org
    else
      dir.s = ∇a_org * dir.x - rhs.primal_r
      dir.y = ( rhs.comp_r - dir.s .* y_org ) ./ s_org # (mu_target - s_cur .* y_cur)
    end

    update_kkt_error!(kkt_solver, Inf)
end
