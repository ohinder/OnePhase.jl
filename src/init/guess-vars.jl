function mehortra_guess( nlp, pars, timer, x::Vector, a::Vector, J, g::Vector )
    s = deepcopy(a);
    m = length(s)

    if true
      println("estimating intial y")
      start_advanced_timer(timer, "INIT/estimate_y_tilde")
      y = estimate_y_tilde( J, g, pars )
      if isbad(y)
        warning("y is bad!!!!!")
      end
      pause_advanced_timer(timer, "INIT/estimate_y_tilde")
    else
      y = ones(m)
    end



    if true
      start_advanced_timer(timer, "INIT/mehortra_guarding")

      ais = cons_indicies(nlp)
      bis = bound_indicies(nlp)
      @show ais, bis

      if isbad(s)
        throw(Eval_NaN_error(getbad(s),x,"s"))
      end

      if isbad(y)
        throw(Eval_NaN_error(getbad(y),x,"y"))
      end

      threshold = 1e-8
      s_new, y = mehortra_guarding( deepcopy(s), deepcopy(y), threshold )
      if isbad(s_new)
        throw(Eval_NaN_error(getbad(s_new),x,"s"))
      end

      if pars.init.start_satisfying_bounds
        s[ais] = s_new[ais]
      else
        s = s_new
      end
    end

    y = min.(ones(m) * 1e3, max.(y, 0.1 * ones(m)))

    #@show y, J

    if isbad(s)
      throw(Eval_NaN_error(getbad(s),x,"s"))
    end

    ifree = _i_not_fixed(nlp.nlp)
    uvar = nlp.nlp.meta.uvar[ifree]
    lvar =  nlp.nlp.meta.lvar[ifree]

    @assert(length(bis) == sum(uvar .< Inf) + sum(lvar .> -Inf))
    #lb()

    for i in bound_indicies(nlp)
        if (s[i] <= 0.0)
            @show i, s[i], a[i] #lvar[i], uvar[i]
        end
        if pars.init.start_satisfying_bounds
          @assert(s[i] == a[i])
        end
    end
    pause_advanced_timer(timer, "INIT/mehortra_guarding")

    return s, y
end


function estimate_y_tilde( J::SparseMatrixCSC{Float64,Int64}, g::Array{Float64,1}, pars::Class_parameters )
    try
      println("estimate y_tilde ...")
      @show densest_col(J)
      @show densest_row(J)
      #@time H = Symmetric(J * J' + norm(J,Inf) * 1e-4 * speye( size(J,1) ))
      n = size(J,2); m = size(J,1);

      lambda = 1e-4
      @show lambda
      if false
        @time H = [speye(n) -J'; J lambda * speye(m)]
        @time M = lufact( H );
        rhs = [-g; zeros(m)];
        sol = M \ rhs
        y = sol[(n+1):end]
      else # cholesky factor
        scaling = 1.0 #norm(J,Inf)^2
        @time H = lambda * speye(n) + J' * J / scaling;
        @time M = cholfact( H );
        dx = scaling * (M \ -g)
        y = -J * dx
      end

      println("linear system solved")
      return y
    catch (e)
      println("error in estimate_y_tilde")
      @show e
      return ones(size(J,1))
    end
end

function mehortra_guarding( s_tilde::Array{Float64,1}, y_tilde::Array{Float64,1}, threshold::Float64 )
    δ_s = max(-2.0 * minimum(s_tilde), threshold)
    δ_y = max(-2.0 * minimum(y_tilde), threshold)

    δ_s_tilde = δ_s + 0.5 * dot(s_tilde + δ_s, y_tilde + δ_y) / sum(y_tilde + δ_y)
    δ_y_tilde = δ_y + 0.5 * dot(s_tilde + δ_s, y_tilde + δ_y) / sum(s_tilde + δ_s)

    s = s_tilde + δ_s_tilde
    y = y_tilde + δ_y_tilde

    return s, y
end
#=
kkt_solver = pick_KKT_solver(pars);
form_system!(kkt_solver, iter, timer)
fact_succeed, inertia_num_fac, new_delta = ipopt_strategy!(iter, kkt_solver, pars, timer)
kkt_associate_rhs!(kkt_solver, iter, Reduct_affine(), timer)
compute_direction!(kkt_solver, timer)
=#
function change_mu!(iter::Class_iterate,new_mu::Float64, pars::Class_parameters)
    iter.point.mu = new_mu
    #iter.primal_residual_intial = (get_cons(iter) - iter.point.s) / new_mu
    #iter.primal_residual_intial = (iter.point.s - get_cons(iter)) / new_mu

    y_c = new_mu ./ iter.point.s
    buffer = 2.0
    y = iter.point.y
    y = min.( y_c / (pars.ls.comp_feas * buffer), max.(y, pars.ls.comp_feas * y_c * buffer)) # project onto complementarity constraints

    iter.point.y = y
end

#=
GLOBAL_COUNT = 0
function KNITRO_guess(iter::Class_iterate, dir::Class_point, pars::Class_parameters)
    global GLOBAL_COUNT
    GLOBAL_COUNT += 1
    if GLOBAL_COUNT == 100
      y = iter.point.y
      s = iter.point.s
      y_temp = y + dir.y
      s_temp = s + dir.s
      y_est = y_temp - 2.0 * min(0.0, minimum(y))
      s_est = s_temp - 2.0 * min(0.0, minimum(s))

      mu = iter.point.mu
      mu_est = norm(y,Inf) * norm(get_primal_res(iter),Inf)
      #0.5 * (iter)+ mean(-y_guess .* get_primal_res(iter))
      #mu_est = norm(dir.x,1) #min(mu * 2.0, max(mu/2.0, mu_est)) #+ get_cons(iter) .* iter.point.y
      suggested_mu = min(100.0 * mu,max(mu_est,0.01 * mu)) #
      #@show suggested_mu
      change_mu!(iter,suggested_mu,pars)
    end
end=#
