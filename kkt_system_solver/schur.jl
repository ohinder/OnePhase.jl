type Schur_KKT_solver <: abstract_KKT_system_solver
    ls_solver::abstract_linear_system_solver
    M::SparseMatrixCSC{Float64,Int64}
    factor_it::Class_iterate
    delta::Float64
    rhs::System_rhs
    dir::Class_point
    kkt_err_norm::Class_kkt_error
    rhs_norm::Float64

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

    scaling = zeros(dim(iter))
    for i = 1:dim(iter)
        scaling[i] = norm(∇a[:,i], Inf)
    end

    #@show size(∇a)
    #@show maximum(scaling), minimum(scaling)
    #tic()

    Σ = spdiagm(iter.point.y ./ iter.point.s) #Base.SparseArrays.CHOLMOD.Sparse(spdiagm(it.point.y ./ it.point.s))
    #@show maximum(Σ)
    #sqrt_Σ = sqrt(Σ)
     #Base.SparseArrays.CHOLMOD.Sparse();
    #∇a_hat = sqrt_Σ * ∇a
    M = (∇a' * Σ * ∇a) + eval_lag_hess(iter)
    pause_advanced_timer("SCHUR/form_system");

    start_advanced_timer("SCHUR/allocate");
    kkt_solver.M = M;
    kkt_solver.factor_it = iter;
    pause_advanced_timer("SCHUR/allocate");

end

function factor!(kkt_solver::Schur_KKT_solver, shift::Float64)
    return ls_factor!(kkt_solver.ls_solver, kkt_solver.M + speye(dim(kkt_solver.factor_it)) * shift, dim(kkt_solver.factor_it), 0)
end

function compute_direction!(kkt_solver::Schur_KKT_solver)
    factor_it = kkt_solver.factor_it
    ∇a_org = eval_jac(factor_it);
    y_org = get_y(factor_it);
    s_org = get_s(factor_it);

    rhs = kkt_solver.rhs

    #r1 + ∇a_org' * (( r3 + (r2 .* y_org) ) ./ s_org)
    schur_rhs = rhs.dual_r + ∇a_org' * (( rhs.comp_r + (rhs.primal_r .* y_org) ) ./ s_org);
    #@show schur_rhs

    dir = kkt_solver.dir;
    dir.x = ls_solve(kkt_solver.ls_solver, schur_rhs)[:];
    res1 = schur_rhs - (kkt_solver.M * dir.x + kkt_solver.delta * dir.x)
    dir.x += ls_solve(kkt_solver.ls_solver, res1)[:];
    res2 = schur_rhs - (kkt_solver.M * dir.x + kkt_solver.delta * dir.x)

    alt_dir_x = schur_rhs / kkt_solver.delta
    #@show norm(res1,Inf)
    #@show norm(res2,Inf)
    #@show norm(kkt_solver.M * alt_dir_x + kkt_solver.delta * alt_dir_x - schur_rhs, Inf)
    #@show norm(, Inf)


    #=res1 = schur_rhs - (kkt_solver.M + kkt_solver.delta) * dir.x
    @show norm(res1, Inf) / norm(schur_rhs,Inf)
    dir.x += ls_solve(kkt_solver.ls_solver, res1)[:];
    res2 = schur_rhs - (kkt_solver.M + kkt_solver.delta) * dir.x
    @show norm(res2, Inf) / norm(schur_rhs,Inf)=#

    dir.s = ∇a_org * dir.x - rhs.primal_r
    dir.y = ( rhs.comp_r - dir.s .* y_org ) ./ s_org # (mu_target - s_cur .* y_cur)
end
