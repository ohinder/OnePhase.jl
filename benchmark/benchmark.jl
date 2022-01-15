include("../src/OnePhase.jl")

using JLD, OnePhase, NLPModels, CUTEst
USE_IPOPT = false
if USE_IPOPT
  using Ipopt
end
include("ipopt.jl")
include("plots.jl")
include("summary.jl")
include("process_results.jl")
include("run_and_store.jl")

function if_mkdir(dir::String)
  if !isdir(dir)
     mkdir(dir)
  end
end;
