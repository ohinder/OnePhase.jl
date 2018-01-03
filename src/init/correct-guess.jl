function correct_guess1( nlp, pars, timer, x, a, J, g, s, y )
      start_advanced_timer(timer, "INIT/centre")

      ais = cons_indicies(nlp)

      y = y / pars.mu_scale
      mu = dot(s,y)[1] / length(y)

      for i = 1:10
        @show size(s), size(y)
        @show mu, norm(s,Inf)

        y_c = mu ./ s
        buffer = 2.0
        y = min.( y_c / (pars.ls.comp_feas * buffer), max.(y, pars.ls.comp_feas * y_c * buffer)) # project onto complementarity constraints

        #Delta_s = norm(g - J' * y,1) / (length(y) + norm(y,1))
        Delta_s = norm(g - J' * y,Inf) / (1.0 + norm(y,Inf))

        primal_feas_avg = norm(a - s,1) / length(s)
        @show primal_feas_avg, Delta_s, norm(y,Inf)
        if Delta_s > primal_feas_avg
          if pars.start_satisfying_bounds
            s[ais] = s[ais] + Delta_s
          else
            s = s + Delta_s
          end

          mu *= 10.0
        else
          break
        end
      end

      pause_advanced_timer(timer, "INIT/centre")

      #Delta_s = norm(g - J' * y,Inf) / (1.0 + 0.01 * norm(y,Inf))
      #@show Delta_s

      start_advanced_timer(timer, "INIT/construct_class")
      init_point = Class_point(x, y, s, mu)
      check_for_nan(init_point)
      init_it = Class_iterate(init_point, nlp, Class_local_info(), timer, pars);
      pause_advanced_timer(timer, "INIT/construct_class")

      return init_it
end


function correct_guess2( nlp, pars, timer, x, a, J, g, y, mu, conWeight )
      start_advanced_timer(timer, "INIT/centre")

      ais = cons_indicies(nlp)
      s = zeros(length(a))

      i = 1
      min_s = deepcopy(s)
      for i = 1:100
        @show mu
        s = a + mu * conWeight
        #@assert()
        y_c = mu ./ s
        buffer = 2.0
        y = min.( y_c / (pars.ls.comp_feas * buffer), max.(y, pars.ls.comp_feas * y_c * buffer)) # project onto complementarity constraints

        min_s = (mu / 20.0) * conWeight #+ conWeight * max(0.0,-minimum(a))

        if all(s .>= min_s) && norm(g - J' * y,1) / (length(s) + norm(y,1)) < mu * pars.init.dual_threshold
          break
        else
          mu *= 2.0
        end
      end

      if i == 100
        @show mu, norm(conWeight), all(s .>= min_s), minimum(a)
        @show norm(g - J' * y,1) / (length(s) + norm(y,1))
        error("init error max it")
      end

      pause_advanced_timer(timer, "INIT/centre")

      start_advanced_timer(timer, "INIT/construct_class")
      init_point = Class_point(x, y, s, mu)
      check_for_nan(init_point)
      init_it = Class_iterate(init_point, nlp, Class_local_info(), timer, pars);
      pause_advanced_timer(timer, "INIT/construct_class")

      return init_it
end


function correct_guess3( nlp, pars, timer, x::Vector, a::Vector, J, g::Vector, y::Vector, mu::Float64, conWeight )
      ais = cons_indicies(nlp)
      bis = bound_indicies(nlp)

      #if norm(g - J' * y,1) >= norm(g,1)
      #  y = 1e-10 * ones(length(y))
      #end

      #n = length(x)
      #tau = norm(g - J' * y,1) / ( n * mu)
      #tau = norm(g - J' * y,1) / ( ncon(nlp) * mu)
      #if tau > 1.0
      #  mu = mu * tau
      #end


      s = a + mu * conWeight
      #s[bis] = a[bis]
      #@show bis, ais
      @assert(all(s[bis] .>= 0.0))
      @assert(all(s[ais] .>= 0.0))
      mu *= pars.init.mu_scale



      #@assert(mu >= 0.0)



      start_advanced_timer(timer, "INIT/construct_class")
      init_point = Class_point(x, y, s, mu)
      check_for_nan(init_point)
      init_it = Class_iterate(init_point, nlp, Class_local_info(), timer, pars);
      pause_advanced_timer(timer, "INIT/construct_class")

      center_dual!(init_it, pars)

      return init_it
end
