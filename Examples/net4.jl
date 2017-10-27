include("../include.jl")
using Ipopt

nlp_raw2 = CUTEstModel("NET4")

function compare_objects(obj1,obj2)
  for n in fieldnames(obj1)
     if(getfield(obj1,n) != getfield(obj2,n))
        println(n, ":")
        println(getfield(obj1,n), " != ", getfield(obj2,n))
     end
  end
end

ORG_STDOUT = STDOUT
file = open("net4.txt", "w")
redirect_stdout(file)

IpoptSolve(nlp_raw2);
one_phase_solve(nlp_raw2);

STDOUT = ORG_STDOUT

close(file)

finalize(nlp_raw2)
