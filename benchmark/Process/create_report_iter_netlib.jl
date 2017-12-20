include("../benchmark.jl")
include("create_report.jl")

overlapping_results = get_CUTEst_results()

its, best, ratios, times = compute_its_etc(overlapping_results,MAX_IT=3000);

using PyPlot
PyPlot.close()
plot_iteration_ratios(its, best, ratios)
savefig("$folder/inf_iter_ratios.pdf")
