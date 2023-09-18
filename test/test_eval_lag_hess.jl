function test_eval_lag_hess_no_constraint()
    model_test = Model()
    set_optimizer(model_test, OnePhase.OnePhaseSolver)
    set_optimizer_attribute(model_test, "output_level", 0)
    @variable(model_test, x[1:2] >= 0.0)
    @constraint(model_test, x[1] + x[2] >= 1.0)
    @NLobjective(model_test, Min, x[1] ^ 3 + x[2] ^ 3)


    nlp = MathOptNLPModel(model_test)
    m = OnePhase.Class_CUTEst(nlp)
    x = [0.5, 0.5]
    y = [0.0, 0.0, 0.0]
    w = 1.0

    @test OnePhase.eval_lag_hess(m, x, y, w) == sparse([1, 2], [1, 2], [3.0, 3.0], 2, 2)
end


function test_eval_lag_hess_one_constraint()
    # model_test = Model(Ipopt.Optimizer)
    model_test = Model()
    set_optimizer(model_test, OnePhase.OnePhaseSolver)
    set_optimizer_attribute(model_test, "output_level", 0)

    @variable(model_test, x[i=1:2] >= 0)
    @NLconstraint(model_test, c, (x[1] + x[2]) ^ 2 >= 1.0)
    @NLobjective(model_test, Min, x[1] ^ 3 + x[2] ^ 3)

    nlp = MathOptNLPModel(model_test)
    m = OnePhase.Class_CUTEst(nlp)
    # m = OnePhase.Class_CUTEst((backend(model_test)).optimizer.model)
    w = 1.0
    x = [0.5, 0.5]
    y = [-0.375, 0, 0]

    @test OnePhase.eval_lag_hess(m, x, y, w) == sparse([1, 2, 2], [1, 1, 2], [3.75, 0.75, 3.75], 2, 2)
end
