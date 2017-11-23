function projection_onto_bounds_ipopt_style( nlp::Class_CUTEst, pars::Class_parameters, x::Array{Float64,1} )
    println("projecting onto bounds ...")

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

      p_L = min(κ_1 * max(1, abs(lvar[i])), κ_2 * (uvar[i] - lvar[i]) )
      p_U = min(κ_1 * max(1, abs(uvar[i])), κ_2 * (uvar[i] - lvar[i]) )

      if lvar[i] > -Inf
        b_L[i] = lvar[i] + p_L
      else
        b_L[i] = -Inf
      end

      if uvar[i] < Inf
        b_U[i] = uvar[i] - p_U
      else
        b_U[i] = Inf
      end

      x[i] = x[i] #+ rand() * 100.0

      if b_L[i] < b_U[i]
        x[i] = min(max(b_L[i], x[i]), b_U[i])
      else
        @show b_L[i], b_U[i]
        error("this shouldn't happen!")
        #x[i] = (uvar[i] + lvar[i]) / 2.0
      end


      if(x[i] <= lvar[i] || x[i] >= uvar[i])
         @show (b_L[i], b_U[i]), (lvar[i], uvar[i]) , x[i]
         error("Error in projection_onto_bounds")
      end
    end

    return x
end
