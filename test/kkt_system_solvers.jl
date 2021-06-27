#using SparseArrays
using MathProgBase
using NLPModels, NLPModelsJuMP

function test_compute_indicies()
    @testset "test_compute_indicies" begin
        J = SparseArrays.sparse([[1.0 1.0 1.0];
                    [-1.0 -1.0 -1.0];
                    [1.0 0.0 0.0];
                    [0.0 4.0 3.0];
                    [0.0 2.0 1.5];
                    [3.0 2.0 1.0];
                    [-2.0 -2.0 -2.0]])

        J_T = SparseArrays.sparse(J')

        sorted_cols = OnePhase.sorted_col_list(J_T)
        @test sorted_cols == [4, 5, 3, 6, 1, 2, 7]
        break_points = OnePhase.compute_breakpoints(J_T,sorted_cols)
        @test break_points == [1, 3, 4, 5]

        u = 1.0 ./ collect(1:size(J,1))
        no_para_indicies, para_row_info = OnePhase.compute_indicies(J,u)
        @test no_para_indicies == [1, 3, 4, 6]
        for i = 1:length(no_para_indicies)
            @test para_row_info[i].first == no_para_indicies[i]
        end
        @test para_row_info[1].ls[2].ind == 2
        @test para_row_info[1].ls[2].ratio == -1.0
        @test para_row_info[1].ls[3].ind == 7
        @test para_row_info[1].ls[3].ratio == -2.0
        @test para_row_info[3].ls[2].ind == 5
        @test para_row_info[3].ls[2].ratio == 0.5

        # rows with no duplicates
        @test para_row_info[2].u == u[3]
        @test para_row_info[2].ls[1].g == 1.0
        @test para_row_info[4].u == u[6]
        @test para_row_info[4].ls[1].g == 1.0
    end
end

function test_compute_indicies_speed()
    J = sprand(5000,5000,0.01)
    u = 1.0 ./ collect(1:size(J,1))
    tic()
    no_para_indicies, para_row_info = OnePhase.compute_indicies(J,u)
    toc()
end
#test_compute_indicies_speed()
########################################
# TODO
# test algorithm with various options
# e.g., KKT_solvers
# - schur
# - symmetric
#
#
#####

function test_kkt_solver(jump_model,pars)
    #JuMP.build(jump_model)
    #MOI.Utilities.attach_optimizer(jump_model)
    nlp_raw = OnePhase.MathOptNLPModel(jump_model)
    ##println("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH", nlp_raw)
    #MathProgBase.initialize(eval, [:Grad, :Jac, :Hess, :HessVec, :ExprGraph])
    #onePhaseMod = MOI.RawSolver()
    #onePhaseMod = MOI.get(jump_model, MOI.RawSolver())
    #println("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", onePhaseMod)
    #nlp_raw = OnePhase.MathProgNLPModel(jump_model.internalModel)
    ##println("1 ", bridge_constraints(jump_model))
    ##println("1 ", typeof(jump_model))
    #println("1 ", typeof(onePhaseMod))
    #println("1 ", onePhaseMod.eval == nothing)
    #println("2 ", get_optimizer_attribute(jump_model, "inner"))
    #nlp_raw = OnePhase.MathProgNLPModel(JuMP.internalmodel(jump_model))
    #nlp_raw = OnePhase.MathProgNLPModel(onePhaseMod)
    #nlp_raw = onePhaseMod
    #nlp = OnePhase.Class_CUTEst(nlp_raw)
    nlp = OnePhase.Class_CUTEst(nlp_raw)
    ##println("2HOHOHOHOHOHOHOHHOHOHOH")
    timer = OnePhase.class_advanced_timer()
    OnePhase.start_advanced_timer(timer)

    iter = OnePhase.mehrotra_init(nlp, pars, timer);
    OnePhase.update!(iter, timer, pars) # is this necessary ????

    @show iter.point.x
    @show iter.point.y
    @show iter.point.mu
    @show iter.point.s

    kkt_solver = OnePhase.pick_KKT_solver(pars);
    #println("+++++++++++++++++++++kkt_solver: ", kkt_solver)
    OnePhase.initialize!(kkt_solver, iter)
    OnePhase.form_system!(kkt_solver, iter, timer)
    OnePhase.factor!(kkt_solver, 1e-8, timer)
    reduct_factors = OnePhase.Reduct_affine()
    OnePhase.kkt_associate_rhs!(kkt_solver, iter, reduct_factors, timer)
    OnePhase.compute_direction!(kkt_solver, timer)

    #@show kkt_solver.dir.x
    #@show kkt_solver.dir.y
    #@show kkt_solver.dir.s

    OnePhase.pause_advanced_timer(timer)
    return kkt_solver.dir
end

function test_kkt_solvers(jump_model)
    pars = OnePhase.Class_parameters()
    pars.output_level = 0

#Original Code
    dir_schur = test_kkt_solver(jump_model,pars)
<<<<<<< HEAD
    #@test_broken begin
    begin
=======
#=
    @test_broken begin
>>>>>>> 3ec239fb3ac9f2d2ae0be2c3544248c566977de2
        pars.kkt.kkt_solver_type=:symmetric
        pars.kkt.linear_solver_type=:HSL

        dir_sym = test_kkt_solver(jump_model,pars)

        @test LinearAlgebra.norm(dir_schur.x - dir_sym.x,2)<1e-6
        @test LinearAlgebra.norm(dir_schur.y - dir_sym.y,2)<1e-6
        @test LinearAlgebra.norm(dir_schur.s - dir_sym.s,2)<1e-6
    end
=#
#Updated Code

    begin
        pars.kkt.kkt_solver_type=:symmetric
        pars.kkt.linear_solver_type=:julia

        dir_sym = test_kkt_solver(jump_model,pars)
        println("----------------------------dir_schur.x: ", dir_schur.x)
        println("------------------------------dir_sym.x: ", dir_sym.x)
        println("----------------------------dir_schur.y: ", dir_schur.y)
        println("------------------------------dir_sym.y: ", dir_sym.y)
        println("----------------------------dir_schur.s: ", dir_schur.s)
        println("------------------------------dir_sym.s: ", dir_sym.s)

<<<<<<< HEAD
    #@test_broken begin
    begin
=======
        @test LinearAlgebra.norm(dir_schur.x - dir_sym.x,2)<1e-6
        @test LinearAlgebra.norm(dir_schur.y - dir_sym.y,2)<1e-6
        @test LinearAlgebra.norm(dir_schur.s - dir_sym.s,2)<1e-6
    end

#Original Code
#=
    @test_broken begin
>>>>>>> 3ec239fb3ac9f2d2ae0be2c3544248c566977de2
        pars.kkt.kkt_solver_type=:clever_symmetric
        pars.kkt.linear_solver_type=:HSL

        dir_clever_sym =test_kkt_solver(jump_model,pars)

        @test LinearAlgebra.norm(dir_clever_sym.x - dir_sym.x,2)<1e-6
        @test LinearAlgebra.norm(dir_clever_sym.y - dir_sym.y,2)<1e-6
        @test LinearAlgebra.norm(dir_clever_sym.s - dir_sym.s,2)<1e-6
    end
=#
#Updated Code

    begin
        pars.kkt.kkt_solver_type=:clever_symmetric
        pars.kkt.linear_solver_type=:julia

        dir_clever_sym =test_kkt_solver(jump_model,pars)
        #println("----------------------------dir_schur.x: ", dir_schur.x)
        println("-----------------------dir_clever_sym.x: ", dir_clever_sym.x)
        #println("----------------------------dir_schur.y: ", dir_schur.y)
        println("-----------------------dir_clever_sym.y: ", dir_clever_sym.y)
        #println("----------------------------dir_schur.s: ", dir_schur.s)
        println("-----------------------dir_clever_sym.s: ", dir_clever_sym.s)
        @test LinearAlgebra.norm(dir_schur.x - dir_clever_sym.x,2)<1e-6
        @test LinearAlgebra.norm(dir_clever_sym.y - dir_schur.y,2)<1e-6
        @test LinearAlgebra.norm(dir_clever_sym.s - dir_schur.s,2)<1e-6
    end

end

function test_kkt_solvers()
    @testset "test_kkt_solvers" begin
	jump_model = toy_lp0()
        test_kkt_solvers(jump_model)
#=
        jump_model = toy_lp1()
        test_kkt_solvers(jump_model)

        jump_model = toy_lp2()
        test_kkt_solvers(jump_model)

        jump_model = toy_lp3()
        test_kkt_solvers(jump_model)

        jump_model = toy_lp4()
        test_kkt_solvers(jump_model)

        jump_model = toy_lp5()
        test_kkt_solvers(jump_model)

        jump_model = toy_lp6()
        test_kkt_solvers(jump_model)

        jump_model = toy_lp7()
        test_kkt_solvers(jump_model)

        jump_model = toy_lp8()
        test_kkt_solvers(jump_model)
=#
    end
end

function test_compare_columns()
    @testset "test_compare_columns" begin
        A = SparseArrays.sparse(SparseArrays.sparse([[0 10 0.0]; [1 0 1]; [1 0 0]; [2 0 0];])')
        @test OnePhase.compare_columns(A,1,2) === true
        @test OnePhase.compare_columns(A,2,1) === false
        @test OnePhase.compare_columns(A,2,3) === false
        @test OnePhase.compare_columns(A,3,2) === true
        @test OnePhase.compare_columns(A,3,1) === false
        @test OnePhase.compare_columns(A,1,3) === true
    end
end
