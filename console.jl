using advanced_timer, JLD

MUMPS = false
include("parameters.jl")
include("utils/utils.jl")
include("linear_system_solvers/linear_system_solvers.jl")
include("kkt_system_solver/include.jl")
include("line_search.jl")
include("IPM_tools.jl")
include("toy_examples.jl")
include("one_phase.jl")
include("init.jl")
include("delta_strategy.jl")

# rsync -a --stats /Users/Oliver/Google Drive/Stanford/Research/one-phase-2.0 ohinder@sherlock.stanford.edu:one-phase-2.0
# move inside algorithm/to parameters

my_par = Class_parameters()

if false
reset_advanced_timer()
start_advanced_timer()
#iter, status, hist = run_netlib_lp("AGG",my_par);
#iter, status, hist = run_netlib_lp("PILOT87",my_par);
iter, status, hist = run_netlib_lp("BANDM",my_par);

pause_advanced_timer()
print_timer_stats()
end

problems = CUTEst.select(max_var=100, max_var=1000, min_con=100, max_con=3000)

problems = ["DISCS"]
it_counts = Dict()
for problem_name in problems
      ORG_STDOUT = STDOUT
      file = open("results/log_$(problem_name).txt", "w")
      redirect_stdout(file)
      
      try
        nlp_raw = CUTEstModel(problem_name)
        nlp = Class_CUTEst(nlp_raw)


        intial_it = initial_point(nlp, nlp_raw.meta.x0)
        #intial_it.point.mu = 0.0

        reset_advanced_timer()
        start_advanced_timer()
        iter, status, history, t = one_phase_IPM(intial_it, my_par);
        pause_advanced_timer()
        print_timer_stats()

        save("results/$(problem_name).jld","history",history)

        finalize(nlp_raw)
        close(file)

        it_counts[problem_name] = t;
    catch(e)
      it_counts[problem_name] = NaN;
    end

    redirect_stdout(ORG_STDOUT)
end

save("results/it_counts.jld","it_counts",it_counts)

# INFEASIBLE PROBLEMS
#nlp_raw = CUTEstModel("JUNKTURN")
#nlp_raw = CUTEstModel("DRCAVTY3") # seems to be feasible, IPOPT struggles
#nlp_raw = CUTEstModel("MODEL")
#nlp_raw = CUTEstModel("WOODSNE")
#nlp_raw = CUTEstModel("DISCS") # 2030 ITS, IPOPT THINKS THIS PROBLEM IS INFEASIBLE PROBLEM

if false
using Ipopt
mp = NLPModels.NLPtoMPB(nlp_raw, IpoptSolver())
MathProgBase.optimize!(mp)
end

#nlp_raw = CUTEstModel("PT") # 13 ITS
#nlp_raw = CUTEstModel("AGG") # 153 ITS
#nlp_raw = CUTEstModel("KISSING") # 180 ITS
#nlp_raw = CUTEstModel("KISSING2") # 151 ITS
#nlp_raw = CUTEstModel("FLETCHCR")
#nlp_raw = CUTEstModel("GENHUMPS")
#nlp_raw = CUTEstModel("CATENARY") # >> 100. Infinities!
#nlp_raw = CUTEstModel("ZIGZAG") # 22 ITS
#nlp_raw = CUTEstModel("HVYCRASH") # 39 ITS
#nlp_raw = CUTEstModel("TRAINF") # 140 ITS or 81 ITS if AGG starts
#nlp_raw = CUTEstModel("TRAINH") # >> 1000, STRUGGLING, LINEAR SOLVER IS NOT V. GOOD # IPOPT 58
#nlp_raw = CUTEstModel("ARTIF") # 14 ITS
#nlp_raw = CUTEstModel("AVION2") # HARD and poorly conditioned
#nlp_raw = CUTEstModel("AVGASB") # 9 ITS
#nlp_raw = CUTEstModel("A4X12") # HARD and poorly conditioned



#nlp_raw = CUTEstModel("GENHUMPS")
if true
nlp = Class_CUTEst(nlp_raw)
intial_it = initial_point(nlp, nlp_raw.meta.x0)
#intial_it.point.mu = 0.0

reset_advanced_timer()
start_advanced_timer()
iter, status, hist = one_phase_IPM(intial_it, my_par);
pause_advanced_timer()
print_timer_stats()

finalize(nlp_raw)
end
