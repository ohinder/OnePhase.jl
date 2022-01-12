#cd("/home/fah33/advanced_timer")
#Base.GC.enable(false)
#include("/home/fah33/advanced_timer/src/advanced_timer.jl")
#cd("/home/fah33/OnePhaseOffline")

using JuMP, Test
using SparseArrays
using LinearAlgebra
using Statistics
using Printf

include("../src/OnePhase.jl")
#Uncomment the below three line of codes when working with HSL is needed
OnePhase.setUSE_HSL(true)
@show OnePhase.USE_HSL
OnePhase.loadHSL("../src/linear_system_solvers/")

TEST_MOI = false
include("CUTEst.jl")
include("problems.jl")
include("kkt_system_solvers.jl")
include("linear_system_solvers.jl")
include("test_moi_nlp_solver.jl")

function unit_tests()
    test_compare_columns()
    test_compute_indicies()
    test_linear_solvers()
    test_kkt_solvers()
end

#function basic_tests(solver)
function basic_tests(options::Dict{String, Any})

    @testset "rosenbrook" begin
        test_rosenbrook1(options)
        test_rosenbrook2(options)
        test_rosenbrook3(options)
        test_rosenbrook4(options)
    end

    @testset "LP" begin
        test_toy_lp1(options)
        test_toy_lp2(options)
        test_toy_lp3(options)
        test_toy_lp4(options)
        test_toy_lp5(options)
        test_toy_lp6(options)
        test_toy_lp7(options)
    end
    @testset "infeasible" begin
        model = toy_lp_inf1()
	    attachSolverWithAttributesToJuMPModel(model, options)
        optimize!(model)
        status = MOI.get(model, MOI.TerminationStatus())
        @test status == :Infeasible

        model = toy_lp_inf2()
        attachSolverWithAttributesToJuMPModel(model, options)
        optimize!(model)
        status = MOI.get(model, MOI.TerminationStatus())
        @test status == :Infeasible

        model = circle_nc_inf1()
	    attachSolverWithAttributesToJuMPModel(model, options)
        optimize!(model)
        status = MOI.get(model, MOI.TerminationStatus())
        #println("--------------------------------------", status)
        @test status == :Infeasible
    end

    @testset "convex_nlp" begin
        @testset "circle1" begin
            model = circle1()
            attachSolverWithAttributesToJuMPModel(model, options)
            optimize!(model)
            status = MOI.get(model, MOI.TerminationStatus())

            @test status == :Optimal
            check_circle1(model)
        end
        @testset "circle2" begin
            model = circle2()
            attachSolverWithAttributesToJuMPModel(model, options)
            optimize!(model)
            status = MOI.get(model, MOI.TerminationStatus())

            @test status == :Optimal
            check_circle2(model)
        end

        @testset "quad_opt" begin
            model = quad_opt()
	        attachSolverWithAttributesToJuMPModel(model, options)
            optimize!(model)
            status = MOI.get(model, MOI.TerminationStatus())

            if Test.Pass == @test status == :Optimal
                check_quad_opt(model)
            end
        end
    end

    @testset "nonconvex" begin
        model = circle_nc1()
    	attachSolverWithAttributesToJuMPModel(model, options)
        optimize!(model)
        status = MOI.get(model, MOI.TerminationStatus())

        @test status == :Optimal
        check_circle_nc1(model)

        model = circle_nc2()
	    attachSolverWithAttributesToJuMPModel(model, options)
        optimize!(model)
        status = MOI.get(model, MOI.TerminationStatus())

        @test status == :Optimal
        check_circle_nc2(model)
    end

    @testset "unbounded_opt_val" begin
        model = lp_unbd()
        attachSolverWithAttributesToJuMPModel(model, options)
        optimize!(model)
        status = MOI.get(model, MOI.TerminationStatus())
        @test :Unbounded == status

        model = circle_nc_unbd()
	    attachSolverWithAttributesToJuMPModel(model, options)
        set_optimizer_attribute(model, "output_level", 0)
        optimize!(model)
	    status = MOI.get(model, MOI.TerminationStatus())
        @test status == :Unbounded

        model = quad_unbd()
	    attachSolverWithAttributesToJuMPModel(model, options)
        optimize!(model)
        status = MOI.get(model, MOI.TerminationStatus())
        @test_broken status == :Unbounded
    end

    @testset "unbounded_feasible_region" begin
        # problems with unbounded feasible region
        test_unbd_feas(options)
    end

    @testset "starting point" begin
        test_starting_point(options,0.5)
        test_starting_point(options,-0.5)
    end

end

function basic_tests()
    #max_it = 100
    max_it = 81
    output_level = 0
    a_norm_penalty = 1e-4
    @testset "basic_tests" begin
        @testset "cholseky linear system solve" begin
            #solver = OnePhase.OnePhaseSolver(max_iter=max_it,
            #solver = OnePhase.OnePhaseSolver(term!max_it=max_it,
            #a_norm_penalty = a_norm_penalty,
            #output_level=output_level,
            #kkt!kkt_solver_type=:schur)
            #basic_tests(solver)
	    options = Dict{String, Any}("term!max_it"=>max_it,
	    "a_norm_penalty"=>a_norm_penalty,
            "output_level"=>output_level,
            "kkt!kkt_solver_type"=>:schur)
	    basic_tests(options)
        end
if OnePhase.USE_HSL
        println("HSL not working")
        @testset "Ma97 linear system solve" begin
            #solver = OnePhase.OnePhaseSolver(term!max_it=max_it,
            #a_norm_penalty = a_norm_penalty,
            #output_level=output_level,
            #kkt!kkt_solver_type=:symmetric,
            #kkt!linear_solver_type=:HSL)
            #basic_tests(solver)
	    options = Dict{String, Any}("term!max_it"=>max_it,
	    "a_norm_penalty"=>a_norm_penalty,
            "output_level"=>output_level,
            "kkt!kkt_solver_type"=>:symmetric,
            "kkt!linear_solver_type"=>:HSL)
	    basic_tests(options)
        end
	
        @testset "Ma97 linear system solve with clever elimination" begin
            #solver = OnePhase.OnePhaseSolver(term!max_it=max_it,
            #a_norm_penalty = a_norm_penalty,
            #output_level=output_level,
            #kkt!kkt_solver_type=:clever_symmetric,
            #kkt!linear_solver_type=:HSL)
            #basic_tests(solver)
	    options = Dict{String, Any}("term!max_it"=>max_it,
	    "a_norm_penalty"=>a_norm_penalty,
            "output_level"=>output_level,
            "kkt!kkt_solver_type"=>:clever_symmetric,
            "kkt!linear_solver_type"=>:HSL)
	    basic_tests(options)
        end
end
        @testset "LDLT julia linear system solve" begin
            #solver = OnePhase.OnePhaseSolver(term!max_it=max_it,
            #a_norm_penalty = a_norm_penalty,
            #output_level=output_level,
            #kkt!kkt_solver_type=:clever_symmetric,
            #kkt!linear_solver_type=:julia)
            #basic_tests(solver)
	    options = Dict{String, Any}("term!max_it"=>max_it,
	    "a_norm_penalty"=>a_norm_penalty,
            "output_level"=>output_level,
            "kkt!kkt_solver_type"=>:clever_symmetric,
            "kkt!linear_solver_type"=>:julia)
	    basic_tests(options)
        end

    end
end

function moi_nlp_tests()
    @testset "test_moi_nlp_solver" begin
        test_lp1_feasible_MOI()
        test_lp1_optimal_MOI()
        test_nlp1_feasible_MOI()
        test_nlp1_optimal_MOI()
        test_lp1_feasible_JuMP()
        test_lp1_optimal_JuMP()
        test_nlp1_optimal_JuMP()
    end
end

# lets run the tests!
unit_tests()
moi_nlp_tests()
basic_tests()
#cutest_tests()
#executeCUTEST_Models()
#=
x0 = [-1.2; 1.0]
model = Model() # No solver is required
@variable(model, x[i=1:2], start=x0[i])
@NLobjective(model, Min, (x[1] - 1)^2 + 100 * (x[2] - x[1]^2)^2)
@NLconstraint(model, x[1]^2 + x[2] <= 1.0)
result = OnePhase.one_phase_solve(model)
=#