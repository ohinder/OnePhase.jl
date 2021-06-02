########################
##### FEASIBLE LPS #####
########################

function test_lp1_feasible_MOI()
    solver = OnePhase.OnePhaseSolver()
    MOI.set(solver, MOI.RawParameter("output_level"), 0)

    x = MOI.add_variables(solver, 2)
    l = MOI.add_constraint(solver, x[1], MOI.GreaterThan(0.0))
    #u = MOI.add_constraint(solver, x[2], MOI.LessThan(0.0))
    u = MOI.add_constraint(solver, x[2], MOI.GreaterThan(0.0))
    #e = MOI.add_constraint(solver, x[1] + x[2], MOI.EqualTo(1.0))
    #e = MOI.add_constraint(solver, x[1] + x[2], MOI.LessThan(1.0))
    a = [1.0, 1.0]
    e = MOI.add_constraint(
        solver,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(a, x), 0.0),
        MOI.LessThan(1.0),
    )
    c = [-1.0, -100.0]
    #MOI.set(
    #       solver,
    #       MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
    #       MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(c, x), 0.0),
    #   );
    #MOI.set(solver, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    #@test isnan(MOI.get(solver, MOI.SolveTime()))
    #MOI.optimize!(solver, model)
    #@test MOI.get(solver, MOI.SolveTime()) > 0.0
    #model_test = Model(with_optimizer(OnePhase.OnePhaseSolver))
    model_test = Model()
    set_optimizer(model_test, OnePhase.OnePhaseSolver)
    set_optimizer_attribute(model_test, "output_level", 0)
    @variable(model_test, x >= 0.0)
    #@variable(model_test, y <= 0.0)
    @variable(model_test, y >= 0.0)
    #@NLobjective(model_test, Min, -x - 100 * y)
    @constraint(model_test, x + y <= 1.0)
    MOI.optimize!(solver, model_test)
    status = MOI.get(solver, MOI.TerminationStatus())
    @test status == :Optimal
end

function test_lp1_optimal_MOI()
    solver = OnePhase.OnePhaseSolver()
    MOI.set(solver, MOI.RawParameter("output_level"), 0)

    x = MOI.add_variables(solver, 2)
    l = MOI.add_constraint(solver, x[1], MOI.GreaterThan(0.0))
    #u = MOI.add_constraint(solver, x[2], MOI.LessThan(0.0))
    u = MOI.add_constraint(solver, x[2], MOI.GreaterThan(0.0))
    #e = MOI.add_constraint(solver, x[1] + x[2], MOI.EqualTo(1.0))
    #e = MOI.add_constraint(solver, x[1] + x[2], MOI.LessThan(1.0))
    a = [1.0, 1.0]
    e = MOI.add_constraint(
        solver,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(a, x), 0.0),
        MOI.LessThan(1.0),
    )
    c = [-1.0, -100.0]
    MOI.set(
           solver,
           MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
           MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(c, x), 0.0),
       );
    MOI.set(solver, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    #@test isnan(MOI.get(solver, MOI.SolveTime()))
    #MOI.optimize!(solver, model)
    #@test MOI.get(solver, MOI.SolveTime()) > 0.0
    #model_test = Model(with_optimizer(OnePhase.OnePhaseSolver))
    model_test = Model()
    set_optimizer(model_test, OnePhase.OnePhaseSolver)
    set_optimizer_attribute(model_test, "output_level", 0)
    @variable(model_test, x >= 0.0)
    #@variable(model_test, y <= 0.0)
    @variable(model_test, y >= 0.0)
    #@NLobjective(model_test, Min, -x - 100 * y)
    @constraint(model_test, x + y <= 1.0)
    MOI.optimize!(solver, model_test)
    status = MOI.get(solver, MOI.TerminationStatus())
    #println("##################", MOI.get(solver, MOI.ObjectiveValue()))
    @test status == :Optimal
end

function test_nlp1_feasible_MOI()
    solver = OnePhase.OnePhaseSolver()
    MOI.set(solver, MOI.RawParameter("output_level"), 0)

    x = MOI.add_variables(solver, 2)
    l = MOI.add_constraint(solver, x[1], MOI.GreaterThan(0.0))
    u = MOI.add_constraint(solver, x[2], MOI.GreaterThan(0.0))
    a = [1.0, 1.0]
    e = MOI.add_constraint(
        solver,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(a, x), 0.0),
        MOI.GreaterThan(0.1),
    )
    c = [-1.0, -100.0]
    #MOI.set(
    #       solver,
    #       MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
    #       MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(c, x), 0.0),
    #   );
    #MOI.set(solver, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    #@test isnan(MOI.get(solver, MOI.SolveTime()))
    #MOI.optimize!(solver, model)
    #@test MOI.get(solver, MOI.SolveTime()) > 0.0
    #model_test = Model(with_optimizer(OnePhase.OnePhaseSolver))
    model_test = Model()
    set_optimizer(model_test, OnePhase.OnePhaseSolver)
    set_optimizer_attribute(model_test, "output_level", 0)
    @variable(model_test, x >= 0.0)
    @variable(model_test, y >= 0.0)
    #@NLobjective(model_test, Min, (2.0 - x)^2 + 100 * (y - x^2)^2)
    @constraint(model_test, x + y >= 0.1)
    MOI.optimize!(solver, model_test)
    status = MOI.get(solver, MOI.TerminationStatus())
    @test status == :Optimal
end

function test_nlp1_optimal_MOI()
    solver = OnePhase.OnePhaseSolver()
    MOI.set(solver, MOI.RawParameter("output_level"), 0)

    x = MOI.add_variables(solver, 2)
    l = MOI.add_constraint(solver, x[1], MOI.GreaterThan(0.0))
    u = MOI.add_constraint(solver, x[2], MOI.GreaterThan(0.0))
    a = [1.0, 1.0]
    e = MOI.add_constraint(
        solver,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(a, x), 0.0),
        MOI.GreaterThan(0.1),
    )
    c = [-1.0, -100.0]
    #MOI.set(
    #       solver,
    #       MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
    #       MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(x, x), 0.0),
    #   );
    #MOI.set(solver, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    #@test isnan(MOI.get(solver, MOI.SolveTime()))
    #MOI.optimize!(solver, model)
    #@test MOI.get(solver, MOI.SolveTime()) > 0.0
    #model_test = Model(with_optimizer(OnePhase.OnePhaseSolver))
    model_test = Model()
    set_optimizer(model_test, OnePhase.OnePhaseSolver)
    set_optimizer_attribute(model_test, "output_level", 0)
    @variable(model_test, x >= 0.0)
    @variable(model_test, y >= 0.0)
    @NLconstraint(model_test, x^2 + y^2 <= 1)
    @NLobjective(model_test, Min, x ^ 2 + y ^ 2)
    #@objective(model_test, Min, x + y)
    @constraint(model_test, x + y >= 0.1)
    #@constraint(model_test, x >= 0.0)
    MOI.optimize!(solver, model_test)

    #println("########NLOBJECTIVE##########", MOI.get(solver, MOI.ObjectiveValue()))
    status = MOI.get(solver, MOI.TerminationStatus())
    @test status == :Optimal
    #nlp = MathOptNLPModel(model_test)
    #println("##########################", nlp)
    #pars = OnePhase.create_pars_JuMP(solver.options)
    #println("##########################", pars)
    #iter, status, hist, t, err, timer = OnePhase.one_phase_solve(nlp,pars)
end


function test_lp1_feasible_JuMP()
    #model_test = Model(with_optimizer(OnePhase.OnePhaseSolver))
    model_test = Model()
    set_optimizer(model_test, OnePhase.OnePhaseSolver)
    set_optimizer_attribute(model_test, "output_level", 0)
    @variable(model_test, x >= 0.0)
    #@variable(model_test, y <= 0.0)
    @variable(model_test, y >= 0.0)
    #@NLobjective(model_test, Min, -x - 100 * y)
    @constraint(model_test, x + y <= 1.0)
    JuMP.optimize!(model_test)
    status =  MOI.get(model_test, MOI.TerminationStatus())
    @test status == :Optimal
end

function test_lp1_optimal_JuMP()
    model_test = Model()
    set_optimizer(model_test, OnePhase.OnePhaseSolver)
    set_optimizer_attribute(model_test, "output_level", 0)
    @variable(model_test, x >= 0.0)
    @variable(model_test, y >= 0.0)
    @NLobjective(model_test, Min, -x - 100 * y)
    @constraint(model_test, x + y <= 1.0)
    JuMP.optimize!(model_test)
    status =  MOI.get(model_test, MOI.TerminationStatus())
    #println("##################", JuMP.objective_value(model_test))
    @test status == :Optimal
end

function test_nlp1_optimal_JuMP()
    model_test = Model()
    set_optimizer(model_test, OnePhase.OnePhaseSolver)
    set_optimizer_attribute(model_test, "output_level", 0)
    @variable(model_test, x >= 0.0)
    @variable(model_test, y >= 0.0)
    @NLconstraint(model_test, x^2 + y^2 <= 1)
    @NLobjective(model_test, Min, x ^ 2 + y ^ 2)
    #@NLobjective(model_test, Max, x ^ 2 + y ^ 2)
    @constraint(model_test, x + y >= 0.1)
    JuMP.optimize!(model_test)
    status =  MOI.get(model_test, MOI.TerminationStatus())
    #println("########NLOBJECTIVE##########", JuMP.objective_value(model_test))
    @test status == :Optimal
end

