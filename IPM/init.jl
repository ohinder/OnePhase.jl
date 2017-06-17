function init(nlp::Class_CUTEst, pars::Class_parameters, timer::class_advanced_timer)
    if pars.init_style == :mehotra
      return mehortra_least_squares_estimate( nlp, pars, timer )
    else
      if pars.start_satisfying_bounds
          return initial_point_satisfy_bounds(nlp, pars, timer)
      else
        return initial_point_generic(nlp, pars, timer)
      end
    end
end

function centre_dual!(point::Class_point, buffer::Float64, pars::Class_parameters)
  y_c = point.mu ./ point.s

  point.y = min( y_c / (pars.comp_feas * buffer), max(point.y, pars.comp_feas * y_c * buffer))
end


function init(nlp::Class_QP, pars::Class_parameters)
    x = suggested_starting_point(nlp)


    #feas_mu_ratio = pars.mu_primal_ratio

    #J = eval_jac(nlp, x)

    a = eval_a(nlp, x)
    m = length(a)

    delta = 1e-3
    #s = delta * ones( m ) #max(-a - 1.5 * minimum(a), 1.0)
    s = max(-a, ones( m ) ) #max(-a - 1.5 * minimum(a), 1.0)
    y = ones( m ) #1e-2 ./ s #ones(length(a)) / norm(a,Inf)
    g = eval_grad_lag(nlp, x, y)

    p = Class_point(x, y, s, 1.0)
    iter = Class_iterate(p, nlp, Class_local_info(0.0));

    kkt_solver = pick_KKT_solver(pars);
    initialize!(kkt_solver, iter)
    form_system!(kkt_solver, iter)
    #update_delta!(kkt_solver, 0.0, 0.0)
    factor!(kkt_solver, sqrt(norm(g,1) / length(g)) )
    reduct_factors = Reduct_affine()
    #filter = Clas
    #ipopt_strategy!(iter, reduct_factors, kkt_solver, pars.ls_mode_stable, filter, pars)
    kkt_associate_rhs!(kkt_solver, iter, Reduct_affine())
    compute_direction!(kkt_solver)

    d = kkt_solver.dir

    p.x = p.x + d.x
    #p.y = max(abs(p.y + d.y),1.0)
    #p.s = max(abs(p.s + d.s),1.0)
    #p.s = max(abs(p.s + d.s), minimum(-a))
    threshold = 0.0
    p.s, p.y = mehortra_guarding( p.s + d.s, p.y + d.y, threshold )
    p.mu = mean(p.y .* p.s)

    infeas = norm(eval_a(nlp, x) - p.s,Inf)
    if p.mu > infeas * 1e3
        p.s += p.mu * 1e-3
    end

    if p.mu < infeas * 1e-3
        p.mu = infeas * 1e-3
    end

    centre_dual!(p, 2.0, pars)

    iter.primal_residual_intial = eval_primal_residual(iter)



    #@show p.x
    #@show p.y
    #@show p.s

    return iter #initial_point_generic(nlp, pars)
end

function projection_onto_bounds1( nlp::Class_CUTEst, pars::Class_parameters, x::Array{Float64,1} )
    ifree = _i_not_fixed(nlp.nlp)
    uvar = deepcopy(nlp.nlp.meta.uvar[ifree])
    lvar =  deepcopy(nlp.nlp.meta.lvar[ifree])

    x = deepcopy(x)
    b_L = zeros(length(x))
    b_U = zeros(length(x))

    κ_1 = 1e-2;
    κ_2 = 1e-2;

    threshold_too_close = 1e-8

    for i = 1:length(x)
      if (uvar[i] < lvar[i] + threshold_too_close)
        @show uvar[i] - lvar[i]
        println("ERROR: bounds too close!")
        error("bounds too close!")
      end

      if lvar[i] < Inf
        b_L[i] = lvar[i] + min(κ_1 * max(1, abs(lvar[i])), κ_2 * (uvar[i] - lvar[i]) )
      else
        b_L[i] = -Inf
      end

      if uvar[i] < Inf
        b_U[i] = uvar[i] - min(κ_1 * max(1, abs(uvar[i])), κ_2 * (uvar[i] - lvar[i]) )
      else
        b_U[i] = Inf
      end


      if b_L[i] < b_U[i]
        x[i] = min(max(b_L[i], x[i]), b_U[i])
      else
        x[i] = (uvar[i] + lvar[i]) / 2.0
      end


      if(x[i] <= lvar[i] || x[i] >= uvar[i])
         @show (b_L[i], b_U[i]), (lvar[i], uvar[i]) , x[i]
         error("Error in projection_onto_bounds")
      end
    end

    return x
end

function projection_onto_bounds2( nlp::Class_CUTEst, pars::Class_parameters, x::Array{Float64,1} )
  threshold_too_close = 1e-6

  thres = 1.0;

  ifree = _i_not_fixed(nlp.nlp)
  uvar = deepcopy(nlp.nlp.meta.uvar[ifree])
  lvar =  deepcopy(nlp.nlp.meta.lvar[ifree])

  x = deepcopy(x)

  for i = 1:length(x)
    if (uvar[i] < lvar[i] + threshold_too_close)
      @show uvar[i] - lvar[i]
      println("ERROR: bounds too close!")
      error("bounds too close!")
    end

    if lvar[i] + threshold_too_close > x[i] || x[i] > uvar[i] - threshold_too_close
      if uvar[i] - lvar[i] < thres
        x[i] = (uvar[i] + lvar[i]) / 2.0
      elseif abs(lvar[i]) < abs(uvar[i])
        x[i] = lvar[i] + thres/2.0
      else
        x[i] = uvar[i] - thres/2.0
      end
    end
  end

  return x
end

function isbad(vec)
  for i = 1:length(vec)
    if isbad(vec[i])
      return true
    end
  end

  return false
end

function mehortra_least_squares_estimate( nlp, pars, timer )
    start_advanced_timer(timer, "INIT/x")
    x = suggested_starting_point(nlp)
    pause_advanced_timer(timer, "INIT/x")

    #x += randn(length(x)) * norm(x) / 10.0
    #=
    J_init = eval_jac(nlp, x)
    @assert(!isbad(J_init))
    a_init = eval_a(nlp, x);
    assert(!isbad(a_init))
    =#

    if pars.start_satisfying_bounds
      start_advanced_timer(timer, "INIT/projection_onto_bounds")
      #x = projection_onto_bounds1( nlp, pars, x )
      x = projection_onto_bounds2( nlp, pars, x )
      pause_advanced_timer(timer, "INIT/projection_onto_bounds")
    end

    start_advanced_timer(timer, "INIT/evals")
    a = eval_a(nlp, x);
    s = a;
    m = length(a)
    J = eval_jac(nlp, x)
    g = eval_grad_f(nlp, x)
    @assert(!isbad(J[:]))
    @assert(!isbad(g))
    @assert(!isbad(a))
    pause_advanced_timer(timer, "INIT/evals")

    #@show sum(g)

    #@assert(sum(g) > 0.0)

    start_advanced_timer(timer, "INIT/estimate_y_tilde")
    y = estimate_y_tilde( J, g, pars )
    pause_advanced_timer(timer, "INIT/estimate_y_tilde")

    start_advanced_timer(timer, "INIT/mehortra_guarding")
    threshold = 1e-4 * (max(-minimum(a), 0) + norm(g, Inf))
    @show threshold

    ais = cons_indicies(nlp)
    bis = bound_indicies(nlp)
    @show ais, bis

    #s[ais], y[ais] = mehortra_guarding( deepcopy(s[ais]), deepcopy(y[ais]), threshold )

    s_new, y = mehortra_guarding( deepcopy(s), deepcopy(y), threshold )


    #if (norm(y,Inf) < 1e-1 || norm(y,Inf) > 1e3)
    #  y = ones(m) * 100.0
    #end
    y = min(ones(m) * 1e3, max(y, 0.1 * ones(m)))

    #@show y, J

    if pars.start_satisfying_bounds
      s[ais] = s_new[ais]
    else
      s = s_new
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
        if pars.start_satisfying_bounds
          @assert(s[i] == a[i])
        end
    end
    pause_advanced_timer(timer, "INIT/mehortra_guarding")

    start_advanced_timer(timer, "INIT/centre")
    mu = dot(s,y) / length(y)
    #infeas = norm(a + s,Inf)
    #if mu / infeas > 1e3
    #  mu = infeas * 1e3
    #end

    y_c = mu ./ s
    buffer = 2.0
    y = min( y_c / (pars.comp_feas * buffer), max(y, pars.comp_feas * y_c * buffer)) # project onto complementarity constraints
    pause_advanced_timer(timer, "INIT/centre")

    start_advanced_timer(timer, "INIT/construct_class")
    init_point = Class_point(x, y, s, mu)
    init_it = Class_iterate(init_point, nlp, Class_local_info(0.0), timer);
    pause_advanced_timer(timer, "INIT/construct_class")

    return init_it
end

function estimate_y_tilde( J, g, pars )
    H = Symmetric(J * J' + norm(J,Inf) * 1e-2 * speye( size(J,1) ))
    #@show H, J
    M = cholfact( H );
    rhs = J * g
    #@show rhs
    return M \ rhs
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

function initial_point_satisfy_bounds(nlp::Class_CUTEst, pars::Class_parameters)
    x = suggested_starting_point(nlp)
    x = projection_onto_bounds2( nlp, pars, x )

    feas_mu_ratio = pars.mu_primal_ratio

    m = nbounds_orginal(nlp) + ncons_orginal(nlp)
    init_point = Class_point(x, zeros(m), zeros(m), 0.0)

    a = eval_a(nlp, x);
    J = eval_jac(nlp, x)
    g = eval_grad_lag(nlp, x, zeros(m))

    #mu = (mean(abs(g)) + 1e-3) * infeas * pars.mu_primal_ratio


    #@show mu
    mu = 0.0
    for k = 1:20
        infeas = norm(min(a,0.0),Inf) + mu / pars.mu_primal_ratio
        mu = infeas * pars.mu_primal_ratio

        for i in bound_indicies(nlp)
            @assert(a[i] > 0.0)
            init_point.s[i] = a[i]
            init_point.y[i] = mu / init_point.s[i]
        end

        for i in cons_indicies(nlp)
            init_point.s[i] = max(a[i], 0.0 ) + infeas
            init_point.y[i] = mu / init_point.s[i]
        end

        if true
          dual_norm = norm(eval_grad_lag(nlp, init_point.x, init_point.y),Inf)
          if  dual_norm / max(norm(init_point.y,Inf),1.0) < 10 * mu
            break
          else
            mu = max(5.0 * mu, 1e-8)
          end
        else
          break
        end
    end

    init_point.mu = mu

    init_it = Class_iterate(init_point, nlp, Class_local_info(0.0));

    return init_it
end


function initial_point_generic(nlp::abstract_nlp,pars::Class_parameters, x::Array{Float64,1})
    a = eval_a(nlp, x)
    infeas = norm(min(a,0.0),Inf)
    m = length(a)
    init_point = Class_point(x, zeros(m), zeros(m), infeas)

    feas_mu_ratio = pars.mu_primal_ratio

    g = eval_grad_lag(nlp, x, zeros(m))
    mu = infeas

    init_point.mu = mu
    init_point.s = max(0, a) + infeas
    init_point.y = mu ./ init_point.s;

    init_it = Class_iterate(init_point, nlp, Class_local_info(0.0));

    # make sure
    # dual <= primal = mu

    return init_it
end
