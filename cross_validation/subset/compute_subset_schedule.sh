#!/bin/bash
fold=0

while (($fold < 10))
do
    job_name="subset_${fold}"
    # echo $job_name
    qsub -N ${job_name} -v FOLD=${fold} compute_subset.pbs  
    fold=$((fold + 1))
done
