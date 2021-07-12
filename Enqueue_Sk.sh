#!/bin/bash

date

echo "Enqueuing TurboNLM.NSSI project..."
cd "/home/Oleg.Kharlanov/TurboNLM.NSSI"
qsub -l nodes=1:ppn=8 -N tnlm.nssi run_turbonlm_sk.sh

echo "List of jobs by Oleg.Kharlanov:"
qstat | grep Oleg.Kharlanov

