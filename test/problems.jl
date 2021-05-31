######################
####Utility Method####
######################
function attachSolverWithAttributesToJuMPModel(model:: Model, options::Dict{String, Any})
	set_optimizer(model, OnePhase.OnePhaseSolver)
	for (name, value) in options
        sname = string(name)
	set_optimizer_attribute(model, sname, value)
    end
end

######################
##### ROSENBROOK #####
######################
function rosenbrook1()
    model = Model()
    @variable(model, x)
    @variable(model, y)
    @NLobjective(model, Min, (2.0 - x)^2 + 100 * (y - x^2)^2)
    return model
end

#function test_rosenbrook1(solver)
function test_rosenbrook1(options::Dict{String, Any})
    @testset "test_rosenbrook1" begin
        model = rosenbrook1()
        #setsolver(model,solver)
        ##println("############################################################", solver == nothing)
        #model2 = Model(
        #     optimizer_with_attributes(
        #     OnePhase.OnePhaseSolver, "Presolve" => 0, "OutputFlag" => 1
        #     )
        #)
	#model2 = Model(OnePhase.OnePhaseSolver)
	#set_optimizer_attribute(model2, "Presolve", 0)
	#set_optimizer_attribute(model2, "OutputFlag", 1)
        ##println("############################################################", typeof(solver))
        ##set_optimizer(model, OnePhase.OnePhaseSolver)
        ##set_optimizer(model, solver)
	attachSolverWithAttributesToJuMPModel(model, options)
        optimize!(model)
        status = MOI.get(model, MOI.TerminationStatus())
        @test_broken status == :Optimal
        ##println("############################################################")
    end
end

function rosenbrook2()
    model = Model()
    #model = Model(with_optimizer(OnePhase.OnePhaseSolver))
    @variable(model, x >= 0.0)
    @variable(model, y >= 0.0)
    @NLobjective(model, Min, (2.0 - x)^2 + 100 * (y - x^2)^2)
    @constraint(model, x + y >= 0.1)
    return model
end

#function test_rosenbrook2(solver)
function test_rosenbrook2(options::Dict{String, Any})
    @testset "test_rosenbrook2" begin
        model = rosenbrook2()
        #setsolver(model,solver)
        ##set_optimizer(model, OnePhase.OnePhaseSolver)
        ##set_optimizer(model, solver)
	attachSolverWithAttributesToJuMPModel(model, options)
        optimize!(model)
        status = MOI.get(model, MOI.TerminationStatus())
        @test status == :Optimal
        check_rosenbrook(model)
    end
end

function rosenbrook3()
    model = Model()
    @variable(model, x >= 0.0)
    @variable(model, y >= 0.0)
    @NLobjective(model, Min, (2.0 - x)^2 + 100 * (y - x^2)^2)
    @constraint(model, x^2 + y^2 >= 0.5)
    return model
end

#function test_rosenbrook3(solver)
function test_rosenbrook3(options::Dict{String, Any})
    @testset "test_rosenbrook3" begin
        model = rosenbrook3()
        #setsolver(model,solver)
        ##set_optimizer(model, OnePhase.OnePhaseSolver)
        ##set_optimizer(model, solver)
	attachSolverWithAttributesToJuMPModel(model, options)
        optimize!(model)
        status = MOI.get(model, MOI.TerminationStatus())
        @test status == :Optimal
        check_rosenbrook(model)
    end
end

function rosenbrook4()
    model = Model()
    @variable(model, x >= 0.0)
    @variable(model, y >= 0.0)
    @NLobjective(model, Min, (2.0 - x)^2 + 100 * (y - x^2)^2)
    return model
end

#function test_rosenbrook4(solver)
function test_rosenbrook4(options::Dict{String, Any})
    @testset "test_rosenbrook4" begin
        model = rosenbrook4()
        #setsolver(model,solver)
        ##set_optimizer(model, OnePhase.OnePhaseSolver)
        ##set_optimizer(model, solver)
	attachSolverWithAttributesToJuMPModel(model, options)
        optimize!(model)
        status = MOI.get(model, MOI.TerminationStatus())
        @test_broken status == :Optimal
    end
end


function check_rosenbrook(model)
    tol = 1e-3
    @test abs(getvalue(model[:x]) - 2.0) < tol
    @test abs(getvalue(model[:y]) - 4.0) < tol
end

########################
##### FEASIBLE LPS #####
########################

function toy_lp1()
    model = Model()
    @variable(model, x >= 0.0)
    @variable(model, y >= 0.0)
    @NLobjective(model, Min, -x - 100 * y)
    @constraint(model, x + y <= 1.0)
    return model
end

function check_toy_lp1(model)
    tol = 1e-3
    @test abs(getvalue(model[:x]) - 0.0) < tol
    @test abs(getvalue(model[:y]) - 1.0) < tol
end

#function test_toy_lp1(solver)
function test_toy_lp1(options::Dict{String, Any})
    model = toy_lp1()
    #setsolver(model,solver)
    ##set_optimizer(model, OnePhase.OnePhaseSolver)
    ##set_optimizer(model, solver)
    attachSolverWithAttributesToJuMPModel(model, options)
    optimize!(model)
    status = MOI.get(model, MOI.TerminationStatus())
    @test status == :Optimal
    check_toy_lp1(model)
end

function toy_lp2()
    model = Model()
    @variable(model, 0.0 <= x <= 1.0)
    @variable(model, 0.0 <= y <= 1.0)
    @NLobjective(model, Min, -x - 100 * y)
    @constraint(model, x + y <= 2.0)
    return model
end

function check_toy_lp2(model)
    tol = 1e-3
    @test abs(getvalue(model[:x]) - 1.0) < tol
    @test abs(getvalue(model[:y]) - 1.0) < tol
end

#function test_toy_lp2(solver)
function test_toy_lp2(options::Dict{String, Any})
    model = toy_lp2()
    #setsolver(model,solver)
    ##set_optimizer(model, OnePhase.OnePhaseSolver)
    ##set_optimizer(model, solver)
    attachSolverWithAttributesToJuMPModel(model, options)
    optimize!(model)
    status = MOI.get(model, MOI.TerminationStatus())
    @test status == :Optimal
    check_toy_lp2(model)
end

function toy_lp3()
    model = Model()
    @variable(model, 0.0 <= x <= 1.0)
    @variable(model, 0.0 <= y <= 1.0)
    @NLobjective(model, Min, x)
    @constraint(model, 1.0 <= x + y <= 2.0)
    return model
end

function check_toy_lp3(model)
    tol = 1e-3
    @test abs(getvalue(model[:x]) - 0.0) < tol
    @test abs(getvalue(model[:y]) - 1.0) < tol
end

#function test_toy_lp3(solver)
function test_toy_lp3(options::Dict{String, Any})
    model = toy_lp3()
    #setsolver(model,solver)
    ##set_optimizer(model, OnePhase.OnePhaseSolver)
    ##set_optimizer(model, solver)
    attachSolverWithAttributesToJuMPModel(model, options)
    optimize!(model)
    status = MOI.get(model, MOI.TerminationStatus())
    @test status == :Optimal
    check_toy_lp3(model)
end

function check_toy_lp4(model)
    tol = 1e-3
    @test abs(getvalue(model[:x]) - 0.0) < tol
    @test abs(getvalue(model[:y]) - 1.0) < tol
end

function toy_lp4()
    model = Model()
    @variable(model, x, lower_bound=0.0, upper_bound=1.0)
    @variable(model, y, lower_bound=0.0, upper_bound=1.0)
    @NLobjective(model, Min, x)
    @constraint(model, 1.0 <= x + y <= 2.0)
    return model
end

function check_toy_lp4(model)
    tol = 1e-3
    @test abs(getvalue(model[:x]) - 0.0) < tol
    @test abs(getvalue(model[:y]) - 1.0) < tol
end

#function test_toy_lp4(solver)
function test_toy_lp4(options::Dict{String, Any})
    model = toy_lp4()
    #setsolver(model,solver)
    ##set_optimizer(model, OnePhase.OnePhaseSolver)
    ##set_optimizer(model, solver)
    attachSolverWithAttributesToJuMPModel(model, options)
    optimize!(model)
    status = MOI.get(model, MOI.TerminationStatus())
    @test status == :Optimal
    check_toy_lp4(model)
end

function toy_lp5()
    model = Model()
    @variable(model, x, lower_bound=0.0, upper_bound=1.0)
    @variable(model, y, lower_bound=0.0, upper_bound=1.0)
    @NLobjective(model, Min, x)
    @constraint(model, x + y == 1.0)
    @constraint(model, x * 32.5 + y * 32.5 == 32.5)
    @constraint(model, 3.0 * x + 3.0 * y <= 3.0)
    return model
end

#function test_toy_lp5(solver)
function test_toy_lp5(options::Dict{String, Any})
    model = toy_lp5()
    #setsolver(model,solver)
    ##set_optimizer(model, OnePhase.OnePhaseSolver)
    ##set_optimizer(model, solver)
    attachSolverWithAttributesToJuMPModel(model, options)
    optimize!(model)
    status = MOI.get(model, MOI.TerminationStatus())
    @test status == :Optimal
    check_toy_lp4(model)
end

function toy_lp6()
    model = Model()
    @variable(model, x, lower_bound=0.0, upper_bound=1.0)
    @variable(model, y, lower_bound=0.0, upper_bound=1.0)
    @NLobjective(model, Min, x)
    @constraint(model, x + y == 1.0)
    @constraint(model, x * 5.5 + y * 5.5 == 5.5)
    return model
end

#function test_toy_lp6(solver)
function test_toy_lp6(options::Dict{String, Any})
    model = toy_lp6()
    #setsolver(model,solver)
    ##set_optimizer(model, OnePhase.OnePhaseSolver)
    ##set_optimizer(model, solver)
    attachSolverWithAttributesToJuMPModel(model, options)
    optimize!(model)
    status = MOI.get(model, MOI.TerminationStatus())
    @test status == :Optimal
    check_toy_lp4(model)
end

function toy_lp7()
    model = Model()
    @variable(model, x, lower_bound=0.0, upper_bound=1.0)
    @variable(model, y, lower_bound=0.0, upper_bound=1.0)
    @NLobjective(model, Min, x)
    @constraint(model, 2.0 * x + y == 1.0)
    return model
end

#function test_toy_lp7(solver)
function test_toy_lp7(options::Dict{String, Any})
    model = toy_lp7()
    #setsolver(model,solver)
    ##set_optimizer(model, OnePhase.OnePhaseSolver)
    ##set_optimizer(model, solver)
    attachSolverWithAttributesToJuMPModel(model, options)
    optimize!(model)
    status = MOI.get(model, MOI.TerminationStatus())
    @test status == :Optimal
    check_toy_lp4(model)
end

function toy_lp8()
    model = Model()
    @variable(model, x, lower_bound=0.0, upper_bound=1.0)
    @variable(model, y, lower_bound=0.0, upper_bound=1.0)
    @NLobjective(model, Min, x)
    @constraint(model, x + y >= 1.0)
    @constraint(model, x * 5.5 + y * 5.5 <= 5.5)
    return model
end

#function test_toy_lp8(solver)
function test_toy_lp8(options::Dict{String, Any})
    model = toy_lp8()
    #setsolver(model,solver)
    ##set_optimizer(model, OnePhase.OnePhaseSolver)
    ##set_optimizer(model, solver)
	attachSolverWithAttributesToJuMPModel(model, options)
    optimize!(model)
    status = MOI.get(model, MOI.TerminationStatus())
    @test status == :Optimal
    check_toy_lp4(model)
end


##########################
##### INFEASIBLE LPS #####
##########################

function toy_lp_inf1()
    model = Model()
    @variable(model, x >= 0.0)
    @variable(model, y >= 0.0)
    @NLobjective(model, Min, x + 100 * y)
    @constraint(model, x + 2 * y <= -1.0)
    return model
end

function toy_lp_inf2()
    model = Model()
    @variable(model, x >= 0.0)
    @variable(model, y >= 0.0)
    @NLobjective(model, Min, x + 100 * y)
    @constraint(model, x + 2 * y <= 2.0)
    @constraint(model, x + 2 * y >= 4.0)
    return model
end

#######################
##### CONVEX NLPS #####
#######################

function circle1()
    model = Model()
    @variable(model, x >= 0.0)
    @variable(model, y >= 0.0)
    @NLobjective(model, Min, x + 100 * y)
    @NLconstraint(model, x^2 + y^2 <= 1.0)
    @NLconstraint(model, (x-2.0)^2 + y^2 <= 1.0)
    return model
end

function check_circle1(model)
    tol = 1e-3
    @test abs(getvalue(model[:x]) - 1.0) < tol
    @test abs(getvalue(model[:y]) - 0.0) < tol
end

function circle2()
    model = Model()
    @variable(model, x >= 0.0)
    @variable(model, y >= 0.0)
    @NLobjective(model, Min, x^3 + y^3)
    @NLconstraint(model, x^2 + y^2 <= 1.0)
    return model
end

function check_circle2(model)
    tol = 1e-2
    @test abs(getvalue(model[:x]) - 0.0) < tol
    @test abs(getvalue(model[:y]) - 0.0) < tol
end

function quad_opt()
    model = Model()
    @variable(model, x)
    @variable(model, y)
    @NLobjective(model, Min, y)
    @NLconstraint(model, x^2 <= y)
    return model
end


function check_quad_opt(model)
    tol = 1e-2
    @test abs(getvalue(model[:x]) - 0.0) < tol
    @test abs(getvalue(model[:y]) - 0.0) < tol
end

##########################
##### NONCONVEX NLPS #####
##########################
function circle_nc1()
    model = Model()
    @variable(model, x >= 0.0)
    @variable(model, y >= 0.0)
    @NLobjective(model, Min, x + 100 * y)
    @NLconstraint(model, x^2 + y^2 == 1.0)
    @NLconstraint(model, (x-2.0)^2 + y^2 == 1.0)
    return model
end

function check_circle_nc1(model)
    tol = 1e-3
    @test abs(getvalue(model[:x]) - 1.0) < tol
    @test abs(getvalue(model[:y]) - 0.0) < tol
end


function circle_nc2()
    model = Model()
    @variable(model, x, start = 1.0)
    @variable(model, y, start = 1.0)
    @NLobjective(model, Min, x)
    @NLconstraint(model, x^2 + y^2 == 1.0)
    return model
end

function check_circle_nc2(model)
    tol = 1e-3
    @test abs(getvalue(model[:x]) + 1.0) < tol
    @test abs(getvalue(model[:y]) - 0.0) < tol
end

function circle_nc_inf1()
    model = Model()
    @variable(model, x, start = 1.0)
    @variable(model, y, start = 1.0)
    @NLobjective(model, Min, x)
    @NLconstraint(model, x^2 + y^2 == 1.0)
    @NLconstraint(model, x^2 + 2 * y^2 == 4.0)
    return model
end

##############################
##### UNBOUNDED PROBLEMS #####
##############################

function lp_unbd()
    model = Model()
    @variable(model, x >= 0.0)
    @variable(model, y)
    @NLobjective(model, Min, -x)
    @NLconstraint(model, x - y <= 1.0)
    return model
end

function circle_nc_unbd()
    model = Model()
    @variable(model, x)
    @variable(model, y)
    @NLobjective(model, Min, x + 0.1 * y)
    @NLconstraint(model, x^2 + y^2 >= 1.0)
    return model
end

function quad_unbd()
    model = Model()
    @variable(model, x)
    @variable(model, y)
    @NLobjective(model, Min, x)
    @NLconstraint(model, x^2 <= y)
    return model
end

#####################################
##### UNBOUNDED FEASIBLE REGION #####
#####################################

function unbd_feas()
    model = Model()
    @variable(model, x >= 0.0)
    @variable(model, y >= 0.0)
    @variable(model, z >= 0.0)
    @NLobjective(model, Min, y)
    @NLconstraint(model, x^2 <= y)
    @NLconstraint(model, z >= 0.0)
    return model
end

#function test_unbd_feas(solver)
function test_unbd_feas(options::Dict{String, Any})
    println("test_unbd_feas")
    model = unbd_feas()
    #setsolver(model,solver)
    ##set_optimizer(model, OnePhase.OnePhaseSolver)
    ##set_optimizer(model, solver)
	attachSolverWithAttributesToJuMPModel(model, options)
    optimize!(model)
    status = MOI.get(model, MOI.TerminationStatus())
    @test status == :Optimal
    @test getvalue(model[:z]) < 1e5
end

###########################
##### STARTING POINTS #####
###########################

function starting_point_prob(start::Float64)
    model = Model()
    @variable(model, x, start = start)
    @NLobjective(model, Min, -x^2)
    @NLconstraint(model, -1.0 <= x <= 1.0)
    return model
end

#function test_starting_point(solver,starting_point::Float64)
function test_starting_point(options::Dict{String, Any},starting_point::Float64)
    if starting_point == 0.0
        warn("don't select this as a starting point")
    end
    model = starting_point_prob(starting_point)
    #setsolver(model,solver)
    ##set_optimizer(model, OnePhase.OnePhaseSolver)
    ##set_optimizer(model, solver)
    attachSolverWithAttributesToJuMPModel(model, options)
    optimize!(model)
    status = MOI.get(model, MOI.TerminationStatus())
    @test status == :Optimal
    if sign(starting_point) < 0.0
        @test abs(getvalue(model[:x]) - 1.0) < 1e-4
    else
        @test abs(getvalue(model[:x]) + 1.0) < 1e-4
    end
end