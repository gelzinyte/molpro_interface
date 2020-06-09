#!/bin/bash
#$ -pe smp 4 # no of threads to use
#$ -l h_rt=1:00:00 #max time of job; will get cleared after
#$ -q  'orinoco' # the two queues
#$ -S /bin/bash
#$ -N four
#$ -j yes
#$ -cwd
echo "-- New Job --"
export  OMP_NUM_THREADS=${NSLOTS} # tells not to use more threads than specified above.
#export OMP_NUM_THREADS=2
source /home/eg475/programs/miniconda3/etc/profile.d/conda.sh
conda activate general
echo 'running script'
python call_calculator.py 
echo "-- The End--"
