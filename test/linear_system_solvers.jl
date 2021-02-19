function test_julia_unsym(A,b,n,m,inertia,timer)
    # just one function for all julia solvers.
    safe = false
    recycle = false
    solver_julia = OnePhase.linear_solver_JULIA(:unsymmetric, safe, recycle)
    OnePhase.initialize!(solver_julia)
    @test inertia == OnePhase.ls_factor!(solver_julia, A, n, m, timer)
    res1 = zeros(length(b))
    OnePhase.ls_solve!(solver_julia, b, res1, timer)
    res2 = OnePhase.ls_solve(solver_julia, b, timer)
    @test res1 == res2
    return res1
end

function test_julia_sym(A,b,n,m,inertia,timer)
    safe = false
    recycle = false
    solver_julia = OnePhase.linear_solver_JULIA(:symmetric, safe, recycle)
    OnePhase.initialize!(solver_julia)
    @test inertia == OnePhase.ls_factor!(solver_julia, A, n, m, timer)
    res1 = zeros(length(b))
    OnePhase.ls_solve!(solver_julia, b, res1, timer)
    res2 = OnePhase.ls_solve(solver_julia, b, timer)
    @test res1 == res2
    return res1
end

function test_julia_chol(A,b,n,m,inertia,timer)
    safe = false
    recycle = false
    solver_julia = OnePhase.linear_solver_JULIA(:definite, safe, recycle)
    OnePhase.initialize!(solver_julia)
    @test inertia == OnePhase.ls_factor!(solver_julia, A, n, m, timer)
    res1 = zeros(length(b))
    OnePhase.ls_solve!(solver_julia, b, res1, timer)
    res2 = OnePhase.ls_solve(solver_julia, b, timer)
    @test res1 == res2
    return res1
end

function test_ma57(A,b,n,m,inertia,timer)
    safe = false
    recycle = false
    solver_julia = OnePhase.linear_solver_HSL(:symmetric, safe, recycle)
    OnePhase.initialize!(solver_julia)
    @test inertia == OnePhase.ls_factor!(solver_julia, A, n, m, timer)
    res1 = zeros(length(b))
    OnePhase.ls_solve!(solver_julia, b, res1, timer)
    res2 = OnePhase.ls_solve(solver_julia, b, timer)
    @test res1 == res2
    return res1
end


function run_linear_solvers(A,b,n,m,inertia)
    timer = OnePhase.class_advanced_timer()
    OnePhase.start_advanced_timer(timer)

    tol = 1e-9

    #dir_julia_unsym = test_julia_unsym(A,b,n,m,inertia,timer)
    dir_julia_sym = test_julia_sym(A,b,n,m,inertia,timer)
    dir_julia_chol = test_julia_chol(A,b,n,m,inertia,timer)
    @test norm(dir_julia_sym - dir_julia_chol) < tol

    @test_broken dir_ma57 = test_ma57(A,b,n,m,inertia,timer)
    @test_broken norm(dir_ma57 - dir_julia_chol) < tol

    A_2 = (A + A')
    for i = 1:size(A,1)
        A_2[i,i] = A[i,i]
    end
    #display(full(A))
    #display(full(A_2))

    dir_julia_sym_2 = test_julia_sym(A_2,b,n,m,inertia,timer)
    @test norm(dir_julia_sym - dir_julia_sym_2) < tol
    dir_julia_chol_2 = test_julia_chol(A_2,b,n,m,inertia,timer)
    @test norm(dir_julia_chol - dir_julia_chol_2) < tol
    @test_broken dir_ma57_2 = test_ma57(A_2,b,n,m,inertia,timer)
    @test_broken norm(dir_ma57 - dir_ma57_2) < tol

    OnePhase.pause_advanced_timer(timer)
end

function test_linear_solvers()
    @testset "test_linear_solvers" begin
        n = 10
        m = 0
        A = speye(10)
        b = rand(10)
        inertia = 1
        run_linear_solvers(A,b,n,m,inertia)

        # add some off diagonal elements
        n = 10
        m = 0
        A = speye(10)
        A[10,1] = 0.1
        A[9,2] = 0.1
        b = rand(10)
        inertia = 1
        run_linear_solvers(A,b,n,m,inertia)
    end
end
