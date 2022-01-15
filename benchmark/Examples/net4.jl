using Ipopt
include("../include.jl")

nlp_raw2 = CUTEstModel("NET4")

my_pars = Class_parameters()
one_phase_solve(nlp_raw2, my_pars);
finalize(nlp_raw2)


ORG_STDOUT = STDOUT
file = open("net4.txt", "w")
redirect_stdout(file)

IpoptSolve(nlp_raw2);
my_pars = Class_parameters()
one_phase_solve(nlp_raw2, my_pars);

STDOUT = ORG_STDOUT

close(file)

finalize(nlp_raw2)
