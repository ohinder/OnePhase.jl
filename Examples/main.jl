include("../include.jl")
#ssh ?????
#include("../include.jl")


# LARGE dual variables
#nlp_raw2 = CUTEstModel("HVYCRASH")
#nlp_raw2 = CUTEstModel("MSS1")
# INFEASIBLE PROBLEMS
#nlp_raw2 = CUTEstModel("10FOLDTR")
#nlp_raw2 = CUTEstModel("NCVXQP8","-param","N=1000")
#nlp_raw2 = CUTEstModel("JUNKTURN")
#nlp_raw2 = CUTEstModel("CATENARY")
#nlp_raw2 = CUTEstModel("DRCAVTY3") # seems to be feasible, IPOPT struggles
#nlp_raw2 = CUTEstModel("MODEL")
#nlp_raw2 = CUTEstModel("FLOSP2HL")
#nlp_raw2 = CUTEstModel("FLOSP2HH")
#nlp_raw2 = CUTEstModel("WOODSNE")
#nlp_raw2 = CUTEstModel("FLOSP2HM")
#nlp_raw2 = CUTEstModel("CHNRSBNE")
#nlp_raw2 = CUTEstModel("KTMODEL")
#nlp_raw2 = CUTEstModel("CONT6-QQ")
#nlp_raw2 = CUTEstModel("WACHBIEG")
#nlp_raw2 = CUTEstModel("SPANHYD")
#nlp_raw2 = CUTEstModel("HS12")
#nlp_raw2 = CUTEstModel("OSCIPANE")
#nlp_raw2 = CUTEstModel("TRO4X4")
#nlp_raw2 = CUTEstModel("A4X12")
#nlp_raw2 = CUTEstModel("BRAINPC1")
#nlp_raw2 = CUTEstModel("BRAINPC7")
#nlp_raw2 = CUTEstModel("SYNPOP24")
#nlp_raw2 = CUTEstModel("ACOPR300")
#nlp_raw2 = CUTEstModel("SPIN2OP")
#nlp_raw2 = CUTEstModel("AIRPORT")
#nlp_raw2 = CUTEstModel("CONT6-QQ")
#nlp_raw2 = CUTEstModel("HAIFAL")
#nlp_raw2 = CUTEstModel("HELSBY")
#nlp_raw2 = CUTEstModel("EIGENCCO")
#nlp_raw2 = CUTEstModel("QPCBOEI1")
#nlp_raw2 = CUTEstModel("PT") # 13 ITS
#nlp_raw2 = CUTEstModel("COSHFUN") # 153 ITS
#nlp_raw2 = CUTEstModel("KISSING") # 180 ITS
#nlp_raw2 = CUTEstModel("KISSING2") # 151 ITS
#nlp_raw2 = CUTEstModel("FLETCHCR")
#nlp_raw2 = CUTEstModel("GENHUMPS")
#nlp_raw2 = CUTEstModel("ZIGZAG") # 22 ITS
#nlp_raw2 = CUTEstModel("TRAINF") # 140 ITS or 81 ITS if AGG starts
#nlp_raw2 = CUTEstModel("ARTIF") # 14 ITS, IPOPT infeasible
#nlp_raw2 = CUTEstModel("AVGASB") # 9 ITS
#nlp_raw2 = CUTEstModel("HYDROELM")
#nlp_raw2 = CUTEstModel("STEENBRC")
#nlp_raw2 = CUTEstModel("QPCSTAIR")
#nlp_raw2 = CUTEstModel("QPNBOEI2")
#nlp_raw2 = CUTEstModel("READING1")
#nlp_raw2 = CUTEstModel("YORKNET")
#nlp_raw2 = CUTEstModel("CHAIN")
#nlp_raw2 = CUTEstModel("DRUGDISE")

#nlp_raw2 = CUTEstModel("ROCKET")
#nlp_raw2 = CUTEstModel("DISC2")
#nlp_raw2 = CUTEstModel("OET7")
#nlp_raw2 = CUTEstModel("ACOPR57")
#nlp_raw2 = CUTEstModel("NET3")
#nlp_raw2 = CUTEstModel("NET4") #,"-param","TIME=144")
#nlp_raw2 = CUTEstModel("GPP","-param","N=2000")
#nlp_raw2 = CUTEstModel("ANTWERP")
#nlp_raw2 = CUTEstModel("HYDCAR20")
#nlp_raw2 = CUTEstModel("STEENBRE")
#nlp_raw2 = CUTEstModel("STEENBRC")
#nlp_raw2 = CUTEstModel("EXPFITC")
#nlp_raw2 = CUTEstModel("LAUNCH")
#nlp_raw2 = CUTEstModel("TRIMLOSS")

#nlp_raw2 = CUTEstModel("HAIFAM")
#nlp_raw2 = CUTEstModel("ACOPR118")
#nlp_raw2 = CUTEstModel("LAKES")
#nlp_raw2 = CUTEstModel("ELATTAR")
#nlp_raw2 = CUTEstModel("A4X12")

#nlp_raw2 = CUTEstModel("CRESC50")
#nlp_raw2 = CUTEstModel("LHAIFAM")
#nlp_raw2 = CUTEstModel("MPC10")
#nlp_raw2 = CUTEstModel("ELEC")
#nlp_raw2 = CUTEstModel("CRESC50")
#nlp_raw2 = CUTEstModel("TRO4X4")
#nlp_raw2 = CUTEstModel("METHANL8")
#nlp_raw2 = CUTEstModel("GROUPING")
#nlp_raw2 = CUTEstModel("ARWHDNE")
#nlp_raw2 = CUTEstModel("HYDCAR6")
#nlp_raw2 = CUTEstModel("EXPFITA")
#nlp_raw2 = CUTEstModel("EXPFITC")
#nlp_raw2 = CUTEstModel("MINPERM")
#HVYCRASH,
#DISCS, EQC, HIMMELBJ,  PFIT1, PFIT3, SSEBNLN
#nlp_raw2 = CUTEstModel("ACOPR57")
#nlp_raw2 = CUTEstModel("DISCS")
#nlp_raw2 = CUTEstModel("HYDCAR6")
#nlp_raw2 = CUTEstModel("DISC2")
#nlp_raw2 = CUTEstModel("SSEBNLN")
#nlp_raw2 = CUTEstModel("BATCH")
#nlp_raw2 = CUTEstModel("SAWPATH")
#nlp_raw2 = CUTEstModel("DRCAVTY3")
#nlp_raw2 = CUTEstModel("EIGENC2")

#nlp_raw2 = CUTEstModel("MANNE")
#nlp_raw2 = CUTEstModel("SPINOP")
#nlp_raw2 = CUTEstModel("SYNPOP24")
#nlp_raw2 = CUTEstModel("GPP")
#mean(abs(grad(nlp_raw2, nlp_raw2.meta.x0))) #, maximum(abs(jac(nlp_raw2, nlp_raw2.meta.x0)))
#maximum(abs(grad(nlp_raw2, nlp_raw2.meta.x0)))
#mean(abs(jac(nlp_raw2, nlp_raw2.meta.x0)))
#finalize(nlp_raw2)

#nlp_raw2 = CUTEstModel("MANNE")
## HARD PROBLEMS
#nlp_raw2 = CUTEstModel("QPNBOEI1")
#nlp_raw2 = CUTEstModel("STEENBRD")
#nlp_raw2 = CUTEstModel("ACOPR14")
#nlp_raw2 = CUTEstModel("ACOPR118")
#nlp_raw2 = CUTEstModel("ACOPP57")
#nlp_raw2 = CUTEstModel("ACOPP300")
#nlp_raw2 = CUTEstModel("LEAKNET")
#nlp_raw2 = CUTEstModel("TFI1")
#nlp_raw2 = CUTEstModel("TRAINH") # >> 1000, STRUGGLING, LINEAR SOLVER IS NOT V. GOOD # IPOPT 58
#nlp_raw2 = CUTEstModel("AVION2") # HARD and poorly conditioned
#nlp_raw2 = CUTEstModel("A4X12") # HARD and poorly conditioned
#nlp_raw2 = CUTEstModel("CRESC100") # >> 100. Infinities!
#nlp_raw2 = CUTEstModel("CRESC50")
#nlp_raw2 = CUTEstModel("AVION2")
#nlp_raw2 = CUTEstModel("QPNSTAIR")
#nlp_raw2 = CUTEstModel("YORKNET")
#nlp_raw2 = CUTEstModel("TOYSARAH")
#nlp_raw2 = CUTEstModel("AGG")
#problem_list = get_problem_list(100000,9999999999)
# BIG PROBLEMS
#BDRY2
#YATP2SQ
#OSCIGRNE (SEEMS KINDA DUMB)
#nlp_raw2 = CUTEstModel("BA-L52")
#nlp_raw2 = CUTEstModel("RDW2D51F","-param","N=512")
#nlp_raw2 = CUTEstModel("YATP2SQ","-param","N=200")

#nlp_raw2 = CUTEstModel("TWOD","-param","N=2")

function compare_objects(obj1,obj2)
  for n in fieldnames(obj1)
     if(getfield(obj1,n) != getfield(obj2,n))
        println(n, ":")
        println(getfield(obj1,n), " != ", getfield(obj2,n))
     end
  end
end
nlp_raw2 = CUTEstModel("DITTERT")

#IpoptSolve(nlp_raw2);
my_pars = Class_parameters()
one_phase_solve(nlp_raw2,my_pars);

finalize(nlp_raw2)
