#!/bin/bash
fold=0
start=0
total=4836
while (($start < $total))
do
    end=$((start + 200))
    if (($end >= $total)); then
        end=$((total - 1))
    fi
    job_name="dock_${fold}_${start}_${end}"
    #echo $job_name
    qsub -N ${job_name} -v FOLD=${fold},START=${start},END=${end} smina_model_dock.pbs
    #echo "------------------"
    start=$((end + 1))
done
