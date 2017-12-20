using JLD
USE_IPOPT = false
if USE_IPOPT
  using Ipopt
end
include("../src/OnePhase.jl")
include("ipopt.jl")
include("plots.jl")
include("summary.jl")
include("process_results.jl")
include("run_and_store.jl")
