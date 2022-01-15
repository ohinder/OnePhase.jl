include("include.jl")
include("read_lp.jl")

####################################
### PLOT A TOY PROBLEM
####################################
using JuMP, Ipopt, PyPlot
m = Model()

A = [[1.0 1.0 0.0]; [1.0 0.0 -1.0]]
b, c, lbounds, ubounds = [0.0 0.0]', [0.0 0.0 0.0]', [0.0 0.0 0.0]', [1e30 1e30 1e30]'

ls = Dict(1e-3 => "-", 1e-6 => "-.", 1e-9 => "--", 0.0 => ":")
delta_set = [1e-3, 1e-6, 1e-9, 0.0]

ylims = [1e-2, 1e10]

subplot(121)
for delta = delta_set
  max_iter = 100
  hist = run_IPOPT(A, b, c, lbounds, ubounds, max_iter, delta)
  dual_hist = compute_dual_hist(hist)

  semilogy(1:length(dual_hist), dual_hist, color="black", linestyle=ls[delta], label="δ = $delta", basey=10)
end
title("IPOPT")
ax = gca()
ax[:set_ylim](ylims)
xlabel("iteration")
ylabel("maximum dual variable")

legend()


function run_one_phase(jump_model)
    nlp_raw = MathProgNLPModel(jump_model)
    nlp = Class_CUTEst(nlp_raw)
    timer = class_advanced_timer()
    start_advanced_timer(timer)
    start_advanced_timer(timer, "INIT")
    intial_it = init(nlp, my_par, timer)
    pause_advanced_timer(timer, "INIT")
    @assert(is_feasible(intial_it, my_par.ls.comp_feas))
    iter, status, hist, t, err = one_phase_IPM(intial_it, my_par, timer);
    pause_advanced_timer(timer)

    return hist
end

subplot(122)
for delta = delta_set
  jump_model = read_lp_into_JuMP(A, b, c, lbounds, ubounds, delta, true)
  hist = run_one_phase(jump_model)
  dual_hist = compute_dual_hist(hist)

  semilogy(1:length(dual_hist), dual_hist, color="black", linestyle=ls[delta], label="δ = $delta", basey=10)
end
title("Well behaved IPM")
xlabel("iteration")
ax = gca()
ax[:set_ylim](ylims)


#ax[:legend]()[:draggable]()
