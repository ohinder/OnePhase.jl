########################
##### FEASIBLE LPS #####
########################

function test_lp1_feasible_MOI()
    solver = OnePhase.OnePhaseSolver()
    MOI.set(solver, MOI.RawParameter("output_level"), 0)

    x = MOI.add_variables(solver, 2)
    l = MOI.add_constraint(solver, x[1], MOI.GreaterThan(0.0))
    u = MOI.add_constraint(solver, x[2], MOI.GreaterThan(0.0))
    a = [1.0, 1.0]
    e = MOI.add_constraint(
        solver,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(a, x), 0.0),
        MOI.LessThan(1.0),
    )
    c = [-1.0, -100.0]
    model_test = Model()
    set_optimizer(model_test, OnePhase.OnePhaseSolver)
    set_optimizer_attribute(model_test, "output_level", 0)
    @variable(model_test, x >= 0.0)
    @variable(model_test, y >= 0.0)
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
    u = MOI.add_constraint(solver, x[2], MOI.GreaterThan(0.0))
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
    model_test = Model()
    set_optimizer(model_test, OnePhase.OnePhaseSolver)
    set_optimizer_attribute(model_test, "output_level", 0)
    @variable(model_test, x >= 0.0)
    @variable(model_test, y >= 0.0)
    @constraint(model_test, x + y <= 1.0)
    MOI.optimize!(solver, model_test)
    status = MOI.get(solver, MOI.TerminationStatus())
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

    model_test = Model()
    set_optimizer(model_test, OnePhase.OnePhaseSolver)
    set_optimizer_attribute(model_test, "output_level", 0)
    @variable(model_test, x >= 0.0)
    @variable(model_test, y >= 0.0)
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

    model_test = Model()
    set_optimizer(model_test, OnePhase.OnePhaseSolver)
    set_optimizer_attribute(model_test, "output_level", 0)
    @variable(model_test, x >= 0.0)
    @variable(model_test, y >= 0.0)
    @NLconstraint(model_test, x^2 + y^2 <= 1)
    @NLobjective(model_test, Min, x ^ 2 + y ^ 2)
    @constraint(model_test, x + y >= 0.1)
    MOI.optimize!(solver, model_test)

    status = MOI.get(solver, MOI.TerminationStatus())
    @test status == :Optimal
end


function test_lp1_feasible_JuMP()
    #model_test = Model(with_optimizer(OnePhase.OnePhaseSolver))
    model_test = Model()
    set_optimizer(model_test, OnePhase.OnePhaseSolver)
    set_optimizer_attribute(model_test, "output_level", 0)
    @variable(model_test, x >= 0.0)
    @variable(model_test, y >= 0.0)
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
    @constraint(model_test, x + y >= 0.1)
    JuMP.optimize!(model_test)
    status =  MOI.get(model_test, MOI.TerminationStatus())
    @test status == :Optimal
end
