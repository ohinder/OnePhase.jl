function eval_init(nlp::Class_CUTEst, pars::Class_parameters, timer::class_advanced_timer, x)
    start_advanced_timer(timer, "INIT/evals")

    a = eval_a(nlp, x);
    J = eval_jac(nlp, x)
    g = eval_grad_f(nlp, x)

    s = deepcopy(a);
    m = length(a)

    if length(SparseArrays.nonzeros(J)) > 0
      if(isbad(SparseArrays.nonzeros(J)))
          throw(Eval_NaN_error(getbad(SparseArrays.nonzeros(J)),x,"J"))
      end

      if(isbad(a))
          throw(Eval_NaN_error(getbad(a),x,"a"))
      end

      if isbad(g)
        throw(Eval_NaN_error(getbad(g),x,"g"))
      end
    end
    pause_advanced_timer(timer, "INIT/evals")

    return a, J, g
end
