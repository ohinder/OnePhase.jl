#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=CUTEst_one_phase
#################
#a file for job output, you can check job progress, append the job ID with %j to make it unique
#SBATCH --output=../results/screen/log.%j.out
#################
# a file for errors from the job
#SBATCH --error=../results/errors/log.%j.err
#################
#time you think you need; default is 2 hours
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
##SBATCH --qos=long
#SBATCH --time=2880
###2880
#################
#Quality of Service (QOS); think of it as job priority, there is also --qos=long for with a max job length of 7 days, qos normal is 48 hours.
# REMOVE "normal" and set to "long" if you want your job to run longer than 48 hours,
# NOTE- in the hns partition the default max run time is 7 days , so you wont need to include qos
#SBATCH --qos=normal
# We are submitting to the dev partition, there are several on sherlock: normal, gpu, owners, hns, bigmem (jobs requiring >64Gigs RAM)
# The more partitions you can submit to the less time you will wait, you can submit to multiple partitions with -p at once in comma separated format.
# srun -p owners -c 8 --mem=32000 -t 2880 --mail-type=ALL --mail-user=ohinder  --epilog /home/wlaw/killvnc.sh /home/wlaw/sherlock-eznovnc.sh
#SBATCH -p jduchi
#################
#number of nodes you are requesting, the more you ask for the longer you wait
#SBATCH --nodes=1
#################
# --mem is memory per node; default is 4000 MB per CPU, remember to ask for enough mem to match your CPU request, since
# sherlock automatically allocates 8 Gigs of RAM/CPU, if you ask for 8 CPUs you will need 32 Gigs of RAM, so either
# leave --mem commented out or request >= to the RAM needed for your CPU request.
#SBATCH --mem=8000
#################
# Have SLURM send you an email when the job ends or fails, careful, the email could end up in your clutter folder
# Also, if you submit hundreds of jobs at once you will get hundreds of emails.
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
# Remember to change this to your email
#SBATCH --mail-user=ohinder@stanford.edu
#now run normal batch commands
ml load CUTEst/linux-cutest
ml load julia/precompiled/0.5.0
ml load hdf5
# You can use srun if your job is parallel
# otherwise:
julia run_one_phase.jl
julia infeas.jl
