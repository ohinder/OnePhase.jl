kinit ohinder@stanford.edu
ssh ohinder@sherlock.stanford.edu

sdev
ml load CUTEst/linux-cutest
ml load julia/precompiled/0.5.0
ml load hdf5
julia main.jl

# sync to sherlock
rsync -a --stats /Users/Oliver/Google\ Drive/Stanford/Research/one-phase-2.0 ohinder@sherlock.stanford.edu:
one-phase-2.0

# download results
#rsync -a --stats ohinder@sherlock.stanford.edu:one-phase-2.0/results /Users/Oliver/Google\ Drive/Stanford/Research/one-phase-2.0

rsync -a --stats ohinder@sherlock.stanford.edu:one-phase-2.0/results/ /Users/Oliver/Documents/IPM-results


#rsync -a --stats ohinder@sherlock.stanford.edu:one-phase-2.0 /Users/Oliver/Google\ Drive/Stanford/Research/OLD

sbatch run_cutest.sbatch
squeue -u $USER

#cd one-phase-2.0
#julia
sbatch install_ipopt.sbatch

ln -s /Users/Oliver/Documents/IPM-results/ results

Pkg.add("Ipopt")
Pkg.add("JuMP")
Pkg.add("JLD")
Pkg.add("CUTEst")
Pkg.add("Calculus")
Pkg.clone("https://github.com/ohinder/advanced_timer.jl.git")
#Pkg.add("...")
