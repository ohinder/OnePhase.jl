using JuMP, Test
include("../src/OnePhase.jl")
include("problems.jl")
include("kkt_system_solvers.jl")
include("linear_system_solvers.jl")

function unit_tests()
    test_compare_columns()
    test_compute_indicies()
    test_linear_solvers()
    test_kkt_solvers()
end

function basic_tests(solver)
    @testset "rosenbrook" begin
        test_rosenbrook1(solver)
        test_rosenbrook2(solver)
        test_rosenbrook3(solver)
        test_rosenbrook4(solver)
    end

    @testset "LP" begin
        test_toy_lp1(solver)
        test_toy_lp2(solver)
        test_toy_lp3(solver)
        test_toy_lp4(solver)
        test_toy_lp5(solver)
        test_toy_lp6(solver)
        test_toy_lp7(solver)
    end

    @testset "infeasible" begin
        model = toy_lp_inf1()
        setsolver(model,solver)
        @test :Infeasible == solve(model)

        model = toy_lp_inf2()
        setsolver(model,solver)
        @test :Infeasible == solve(model)

        model = circle_nc_inf1()
        setsolver(model,solver)
        @test :Infeasible == solve(model)
    end

    @testset "convex_nlp" begin
        @testset "circle1" begin
            model = circle1()
            setsolver(model,solver)
            @test :Optimal == solve(model)
            check_circle1(model)
        end
        @testset "circle2" begin
            model = circle2()
            setsolver(model,solver)
            @test :Optimal == solve(model)
            check_circle2(model)
        end

        @testset "quad_opt" begin
            model = quad_opt()
            setsolver(model,solver)
            if Base.Test.Pass == @test :Optimal == solve(model)
                check_quad_opt(model)
            end
        end
    end

    @testset "nonconvex" begin
        model = circle_nc1()
        setsolver(model,solver)
        @test :Optimal == solve(model)
        check_circle_nc1(model)

        model = circle_nc2()
        setsolver(model,solver)
        @test :Optimal == solve(model)
        check_circle_nc2(model)
    end

    @testset "unbounded_opt_val" begin
        model = lp_unbd()
        setsolver(model,solver)
        @test :Unbounded == solve(model)

        model = circle_nc_unbd()
        setsolver(model,solver)
        status = solve(model)
        @test status == :Unbounded

        model = quad_unbd()
        setsolver(model,solver)
        status = solve(model)
        @test_broken status == :Unbounded
    end

    @testset "unbounded_feasible_region" begin
        # problems with unbounded feasible region
        test_unbd_feas(solver)
    end

    @testset "starting point" begin
        test_starting_point(solver,0.5)
        test_starting_point(solver,-0.5)
    end
end

function basic_tests()
    max_it = 100
    output_level = 0
    a_norm_penalty = 1e-4
    @testset "basic_tests" begin
        @testset "cholseky linear system solve" begin
            solver = OnePhase.OnePhaseSolver(term!max_it=max_it,
            a_norm_penalty = a_norm_penalty,
            output_level=output_level,
            kkt!kkt_solver_type=:schur)
            basic_tests(solver)
        end

        println("HSL not working")
        #=
        @testset "Ma97 linear system solve" begin
            solver = OnePhase.OnePhaseSolver(term!max_it=max_it,
            a_norm_penalty = a_norm_penalty,
            output_level=output_level,
            kkt!kkt_solver_type=:symmetric,
            kkt!linear_solver_type=:HSL)

            basic_tests(solver)
        end

        @testset "Ma97 linear system solve with clever elimination" begin
            solver = OnePhase.OnePhaseSolver(term!max_it=max_it,
            a_norm_penalty = a_norm_penalty,
            output_level=output_level,
            kkt!kkt_solver_type=:clever_symmetric,
            kkt!linear_solver_type=:HSL)
            basic_tests(solver)
        end
        =#

        @testset "LDLT julia linear system solve" begin
            solver = OnePhase.OnePhaseSolver(term!max_it=max_it,
            a_norm_penalty = a_norm_penalty,
            output_level=output_level,
            kkt!kkt_solver_type=:clever_symmetric,
            kkt!linear_solver_type=:julia)
            basic_tests(solver)
        end
    end
end

# lets run the tests!
unit_tests()
basic_tests()
