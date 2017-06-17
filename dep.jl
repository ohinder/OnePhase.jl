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
rsync -a --stats ohinder@sherlock.stanford.edu:one-phase-2.0/results /Users/Oliver/Google\ Drive/Stanford/Research/one-phase-2.0

#rsync -a --stats ohinder@sherlock.stanford.edu:one-phase-2.0 /Users/Oliver/Google\ Drive/Stanford/Research/OLD

sbatch run_cutest.sbatch
squeue -u $USER

#cd one-phase-2.0
#julia
