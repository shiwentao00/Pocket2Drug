#!/bin/bash
#PBS -q work
#PBS -l nodes=1:ppn=20
#PBS -l walltime=72:00:00
#PBS -A hpc_michal01
#PBS -j oe

cd /work/derick/siamese-monet-project/Pocket2Drug/cross_validation/subset/

singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python compute_subset.py

