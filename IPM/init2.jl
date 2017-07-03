function init(nlp::Class_CUTEst, pars::Class_parameters, timer::class_advanced_timer)
    return mehortra_least_squares_estimate( nlp, pars, timer )
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

function mehortra_least_squares_estimate( nlp, pars, timer )
    start_advanced_timer(timer, "INIT/x")
    x = suggested_starting_point(nlp)
    pause_advanced_timer(timer, "INIT/x")

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

    if length(nonzeros(J)) > 0
      if(isbad(nonzeros(J)))
          throw(Eval_NaN_error(getbad(nonzeros(J)),x,"a"))
      end

      if(isbad(a))
          throw(Eval_NaN_error(getbad(a),x,"a"))
      end

      if isbad(g)
        throw(Eval_NaN_error(getbad(g),x,"g"))
      end
    end
    pause_advanced_timer(timer, "INIT/evals")


    #@show sum(g)

    #@assert(sum(g) > 0.0)

    start_advanced_timer(timer, "INIT/estimate_y_tilde")
    y = estimate_y_tilde( J, g, pars )
    pause_advanced_timer(timer, "INIT/estimate_y_tilde")

    start_advanced_timer(timer, "INIT/mehortra_guarding")
    threshold = 1e-8 #1e-4 * (max(-minimum([a; 0.0]), 0) + norm(g, Inf))
    #@show threshold

    ais = cons_indicies(nlp)
    bis = bound_indicies(nlp)
    @show ais, bis

    #s[ais], y[ais] = mehortra_guarding( deepcopy(s[ais]), deepcopy(y[ais]), threshold )

    Delta_s = norm(g - J' * y,Inf) / (1.0 + 0.01 * norm(y,Inf))
    if pars.start_satisfying_bounds
      s[ais] = s[ais] + Delta_s
    else
      s = s + Delta_s
    end
    if isbad(s)
      throw(Eval_NaN_error(getbad(s),x,"s"))
    end

    if isbad(y)
      throw(Eval_NaN_error(getbad(y),x,"y"))
    end

    s_new, y = mehortra_guarding( deepcopy(s), deepcopy(y), threshold )
    if isbad(s_new)
      throw(Eval_NaN_error(getbad(s_new),x,"s"))
    end

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
    check_for_nan(init_point)
    init_it = Class_iterate(init_point, nlp, Class_local_info(0.0, NaN, :init), timer, pars);
    pause_advanced_timer(timer, "INIT/construct_class")

    return init_it
end

function estimate_y_tilde( J, g, pars )
    try
      H = Symmetric(J * J' + norm(J,Inf) * 1e-4 * speye( size(J,1) ))
      #@show H, J
      M = cholfact( H );
      rhs = J * g
      #@show rhs
      return M \ rhs
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
