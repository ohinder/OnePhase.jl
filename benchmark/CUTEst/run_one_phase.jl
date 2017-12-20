# RUN_LIST as input
include("run_cutest.jl")

problem_list = default_list()

# :dynamic,
#RUN_LIST = ARGS
#RUN_LIST = [:corrections, :step_style, :regularizer, :tol]
RUN_LIST = ["Dec19"]#"LowTol_Dec19",

if "HighTol_Dec19" in RUN_LIST
my_par = Class_parameters()
my_par.term.max_time = 60.0 * 60
my_par.term.max_it = 3000
folder_name = "one_phase/Dec19"
if_mkdir("../results/$folder_name")
run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)
end


if "LowTol_Dec19" in RUN_LIST
my_par = Class_parameters()
my_par.term.max_time = 60.0 * 60
my_par.term.max_it = 3000
my_par.term.tol_opt = 1e-2
my_par.term.tol_inf_1 = 1e-2
my_par.term.tol_inf_2 = 1e-2
folder_name = "one_phase/LowTol_Dec19"
if_mkdir("../results/$folder_name")
run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)
end

if "init" in RUN_LIST
my_par.term.max_time = 10.0 * 60.0
my_par.term.max_it = 1000
if_mkdir("../results/one_phase/init")

my_par.init.mehotra_scaling = false

my_par.init.mu_scale = 0.1
folder_name = "one_phase/init/0.1"
if_mkdir("../results/$folder_name")
run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)

my_par.init.mu_scale = 1.0
folder_name = "one_phase/init/1"
if_mkdir("../results/$folder_name")
run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)

my_par.init.mu_scale = 10.0
folder_name = "one_phase/init/10"
if_mkdir("../results/$folder_name")
run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)

my_par.init.mu_scale = 100.0
folder_name = "one_phase/init/100"
if_mkdir("../results/$folder_name")
run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)

my_par.init.mu_scale = 1000.0
folder_name = "one_phase/init/1000"
if_mkdir("../results/$folder_name")
run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)
end

#problem_list = get_problem_list(100,200)

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

if "penalty" in RUN_LIST
if_mkdir("../results/one_phase/Oct22")
folder_name = "one_phase/Oct22/penalty"
my_par.x_norm_penalty = 1e-8
my_par.a_norm_penalty = 1e-4
my_par.use_prox = true
my_par.use_reg = true
if_mkdir("../results/$folder_name")
run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)

folder_name = "one_phase/Oct22/penalty_no_prox"
my_par.x_norm_penalty = 1e-8
my_par.a_norm_penalty = 1e-4
my_par.use_prox = false
my_par.use_reg = true
if_mkdir("../results/$folder_name")
run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)

folder_name = "one_phase/Oct22/a_penalty_only"
my_par.x_norm_penalty = 1e-30
my_par.a_norm_penalty = 1e-4
my_par.use_prox = true
my_par.use_reg = true
if_mkdir("../results/$folder_name")
run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)

folder_name = "one_phase/Oct22/a_penalty_only_no_prox"
my_par.x_norm_penalty = 1e-30
my_par.a_norm_penalty = 1e-4
my_par.use_prox = false
my_par.use_reg = true
if_mkdir("../results/$folder_name")
run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)
end

if "test" in RUN_LIST
    problem_list = ["PT","AGG","ROBOT"]
    folder_name = "one_phase/test_run"
    if_mkdir("../results/$folder_name")
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)
end

if "plain" in RUN_LIST
    folder_name = "one_phase/plain"
    if_mkdir("../results/$folder_name")
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)
end

if "dynamic" in RUN_LIST
    folder_name = "one_phase/sept_dynamic"
    if_mkdir("../results/$folder_name")
    my_par.adaptive_mu = :paper2
    my_par.primal_bounds_dual_feas = true
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)
    my_par.adaptive_mu = :none
    my_par.primal_bounds_dual_feas = false
end

if "corrections" in RUN_LIST
    #folder_name = "one_phase/sept_1_corrections"
    #if_mkdir("../results/$folder_name")
    #my_par.max_it_corrections = 1
    #run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)

    folder_name = "one_phase/sept_2_corrections"
    if_mkdir("../results/$folder_name")
    my_par.max_it_corrections = 2
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)

    folder_name = "one_phase/sept_3_corrections"
    if_mkdir("../results/$folder_name")
    my_par.max_it_corrections = 3
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)

    folder_name = "one_phase/sept_4_corrections"
    if_mkdir("../results/$folder_name")
    my_par.max_it_corrections = 4
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)
end

if "step_style" in RUN_LIST
    folder_name = "one_phase/sept_max_step_stable"
    if_mkdir("../results/$folder_name")
    my_par.adaptive_mu = :paper
    my_par.max_it_corrections = 3
    ls_mode = :accept_aggressive
    my_par.ls_mode_stable_trust = ls_mode
    my_par.ls_mode_stable_delta_zero = ls_mode
    my_par.ls_mode_stable_correction = ls_mode
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)

    folder_name = "one_phase/sept_log_barrier_stable"
    if_mkdir("../results/$folder_name")
    my_par.adaptive_mu = :paper
    my_par.max_it_corrections = 3
    ls_mode = :accept_stable
    my_par.ls_mode_stable_trust = ls_mode
    my_par.ls_mode_stable_delta_zero = ls_mode
    my_par.ls_mode_stable_correction = ls_mode
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)

    ls_mode = :accept_filter
    my_par.ls_mode_stable_trust = ls_mode
    my_par.ls_mode_stable_delta_zero = ls_mode
    my_par.ls_mode_stable_correction = ls_mode
end

if "regularizer" in RUN_LIST
    folder_name = "one_phase/sept_no_regularizer"
    if_mkdir("../results/$folder_name")
    my_par.use_prox = false
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)
    my_par.use_prox = true
end

if "tol" in RUN_LIST
    folder_name = "one_phase/sept_high_tol"
    if_mkdir("../results/$folder_name")
    my_par.tol = 1e-8
    run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)
end

if "dual_ls_style" in RUN_LIST
  folder_name = "one_phase/dual_ls_0"
  if_mkdir("../results/$folder_name")
  my_par.dual_ls = 0
  run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)


  folder_name = "one_phase/dual_ls_1"
  if_mkdir("../results/$folder_name")
  my_par.dual_ls = 1
  run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)


  folder_name = "one_phase/dual_ls_2"
  if_mkdir("../results/$folder_name")
  my_par.dual_ls = 2
  run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)
end



if "eta_strategy" in RUN_LIST
  #folder_name = "one_phase/eta/mehotra"
  #if_mkdir("../results/$folder_name")
  #my_par.agg_eta = :mehrotra
  #run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)


  folder_name = "one_phase/eta/mehrotra_stb"
  if_mkdir("../results/$folder_name")
  my_par.agg_eta = :mehrotra_stb
  run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)


  folder_name = "one_phase/eta/affine"
  if_mkdir("../results/$folder_name")
  my_par.agg_eta = :affine
  run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)
end
