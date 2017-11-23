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

    y = min(ones(m) * 1e3, max(y, 0.1 * ones(m)))

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
