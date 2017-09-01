include("run_cutest.jl")
TIME_LIMIT = 60.0 * 60
my_par.MAX_TIME = TIME_LIMIT

RUN_LIST = [:plain]

problem_list = get_problem_list(100,10000)
#test_problems(problem_list,1)
if false
problem_list = get_problem_list(6000,7000)
for problem in problem_list
    @show problem
    nlp_raw = CUTEstModel(problem)
    @show nlp_raw.meta.nnzj, nlp_raw.meta.nnzh
    @show nlp_raw.meta.ncon, nlp_raw.meta.nvar
    A = jac(nlp_raw,nlp_raw.meta.x0);
    @show densest_row(A), densest_col(A)

    finalize(nlp_raw)
end
end

if :test in RUN_LIST
    problem_list = ["PT","AGG","ROBOT"]
    folder_name = "one_phase/test_run"
    if_mkdir("../results/$folder_name")
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)
end

if :plain in RUN_LIST
    folder_name = "one_phase/plain"
    if_mkdir("../results/$folder_name")
    my_par.adaptive_mu = :none
    my_par.primal_bounds_dual_feas = false
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)
end

if :dynamic in RUN_LIST
    folder_name = "one_phase/dynamic"
    if_mkdir("../results/$folder_name")
    my_par.adaptive_mu = :paper2
    my_par.primal_bounds_dual_feas = true
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)

    folder_name = "one_phase/nodynamic"
    if_mkdir("../results/$folder_name")
    my_par.adaptive_mu = :none
    my_par.primal_bounds_dual_feas = false
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)
end

if :corrections in RUN_LIST
    folder_name = "one_phase/large_1_corrections"
    if_mkdir("../results/$folder_name")
    my_par.adaptive_mu = :paper
    my_par.max_it_corrections = 1
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)

    folder_name = "one_phase/large_2_corrections"
    if_mkdir("../results/$folder_name")
    my_par.adaptive_mu = :paper
    my_par.max_it_corrections = 2
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)

    folder_name = "one_phase/large_3_corrections"
    if_mkdir("../results/$folder_name")
    my_par.adaptive_mu = :paper
    my_par.max_it_corrections = 3
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)

    folder_name = "one_phase/large_4_corrections"
    if_mkdir("../results/$folder_name")
    my_par.adaptive_mu = :paper
    my_par.max_it_corrections = 4
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)
end

if :step_style in RUN_LIST
    folder_name = "one_phase/large_max_step_stable"
    if_mkdir("../results/$folder_name")
    my_par.adaptive_mu = :paper
    my_par.max_it_corrections = 3
    ls_mode = :accept_aggressive
    my_par.ls_mode_stable_trust = ls_mode
    my_par.ls_mode_stable_delta_zero = ls_mode
    my_par.ls_mode_stable_correction = ls_mode
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)

    folder_name = "one_phase/large_log_barrier_stable"
    if_mkdir("../results/$folder_name")
    my_par.adaptive_mu = :paper
    my_par.max_it_corrections = 3
    ls_mode = :accept_stable
    my_par.ls_mode_stable_trust = ls_mode
    my_par.ls_mode_stable_delta_zero = ls_mode
    my_par.ls_mode_stable_correction = ls_mode
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)
end

if :regularizer in RUN_LIST
    folder_name = "one_phase/regularizer"
    if_mkdir("../results/$folder_name")
    my_par.use_prox = true
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)

    folder_name = "one_phase/no_regularizer"
    if_mkdir("../results/$folder_name")
    my_par.use_prox = false
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)
end

if :tol in RUN_LIST
    folder_name = "one_phase/norm_tol"
    if_mkdir("../results/$folder_name")
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)

    folder_name = "one_phase/high_tol"
    if_mkdir("../results/$folder_name")
    my_par.tol = 1e-8
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)
end
