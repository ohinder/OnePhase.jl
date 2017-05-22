

function initial_point_satisfy_bounds(nlp::Class_CUTEst, pars::Class_parameters)
    x = suggested_starting_point(nlp)

    threshold_too_close = 1e-6
    thres = 1.0;
    feas_mu_ratio = 1e1

    ifree = _i_not_fixed(nlp.nlp)
    uvar = nlp.nlp.meta.uvar[ifree]
    lvar =  nlp.nlp.meta.lvar[ifree]

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

    m = nbounds_orginal(nlp) + ncons_orginal(nlp)
    init_point = Class_point(x, zeros(m), zeros(m), 0.0)

    a = eval_a(nlp, x);
    mu = feas_mu_ratio * max(norm(min(a,0.0),Inf), 1e-3);

    for k = 1:20
        for i in bound_indicies(nlp)
            @assert(a[i] > 0.0)
            init_point.s[i] = a[i]
            init_point.y[i] = mu / init_point.s[i]
        end

        for i in cons_indicies(nlp)
            init_point.s[i] = max(a[i], 0.0 ) + mu / feas_mu_ratio
            init_point.y[i] = mu / init_point.s[i]
        end

        if true
          if norm(eval_grad_lag(nlp, init_point.x, init_point.y),Inf) / max(norm(init_point.y,Inf),1.0) < 10 * mu
            break
          else
            mu *= 5.0
          end
        else
          break
        end
    end

    init_point.mu = mu

    init_it = Class_iterate(init_point, nlp, Class_local_info(0.0));

    return init_it
end


function initial_point_generic(nlp::abstract_nlp,pars::Class_parameters)
    x = suggested_starting_point(nlp)
    a = eval_a(nlp, x)
    infeas = norm(min(a,0.0),Inf)
    init_point = Class_point(x, zeros(length(a)), zeros(length(a)), infeas)

    feas_mu_ratio = 1e1
    mu = feas_mu_ratio * max(infeas, 1e-8)
    for i = 1:20
      init_point.mu = mu
      init_point.s = max(0, a) + mu / feas_mu_ratio
      init_point.y = mu ./ init_point.s;
    end

    init_it = Class_iterate(init_point, nlp, Class_local_info(0.0));

    # make sure
    # dual <= primal = mu

    return init_it
end
