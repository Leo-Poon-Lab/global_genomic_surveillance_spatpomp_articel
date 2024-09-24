#!/bin/bash

## RUN ON HPC

if [ -z "$1" ]
  then
    echo "N_start of jobs not supplied!"
    exit 1
  else
    N_start_jobs=$1
    echo "N_start_jobs: ${N_start_jobs}"
fi

if [ -z "$2" ]
  then
    echo "N_end of jobs not supplied!"
    exit 1
  else
    N_end_jobs=$2
    echo "N_end_jobs: ${N_end_jobs}"
fi

if [ -z "$3" ]
  then
    echo "PARTITION not supplied!"
    exit 1
  else
    PARTITION=$3
    echo "PARTITION: ${PARTITION}"
fi

if [ ${PARTITION} != "amd" ] && [ ${PARTITION} != "condo_amd" ] && [ ${PARTITION} != "intel" ]
  then
    echo "PARTITION not valid!"
    exit 1
fi

if [ -z "$4" ]
  then
    echo "NCORES not supplied!"
    exit 1
  else
    NCORES=$4
    echo "NCORES: ${NCORES}"
fi

if [ -z "$5" ]
  then
    echo "MODEL not supplied!"
    exit 1
  else
    MODEL=$5
    echo "MODEL: ${MODEL}"
fi

if [ ${PARTITION} == "amd" ]
  then
  seq_i=$(seq ${N_end_jobs} -1 ${N_start_jobs})
  # seq_i=$(seq 17 -1 8)
else
  seq_i=$(seq ${N_start_jobs} ${N_end_jobs})
  # seq_i=$(seq 10 11)
fi

for i in ${seq_i}; do
  job_script="job_$i.sh"
  echo "#!/bin/bash" > $job_script
  echo "" >> $job_script
  echo "## NECESSARY JOB SPECIFICATIONS" >> $job_script
  echo "#SBATCH --job-name=Simulate_${MODEL}_${i}" >> $job_script
  echo "#SBATCH --partition=${PARTITION}" >> $job_script
  echo "#SBATCH --qos=normal" >> $job_script
  echo "#SBATCH --time=24:00:00" >> $job_script
  echo "#SBATCH --nodes=1" >> $job_script
  echo "#SBATCH --ntasks=1" >> $job_script
  echo "#SBATCH --cpus-per-task=${NCORES}" >> $job_script
  # echo "#SBATCH --mem=0" >> $job_script
  echo "#SBATCH --mem-per-cpu=3900M" >> $job_script
  echo "#SBATCH --output=%x_%j.out" >> $job_script
  # echo "#SBATCH --error=%x_%j.err" >> $job_script
  echo "" >> $job_script
  echo "## OPTIONAL JOB SPECIFICATIONS" >> $job_script
  echo "#SBATCH --mail-user=guhaogao@hku.hk" >> $job_script
  echo "#SBATCH --mail-type=ALL" >> $job_script
  echo "" >> $job_script
  echo "module load gcc/10.2" >> $job_script
  echo "module load cmake/3.19.7" >> $job_script
  echo "module load libjpeg-turbo/2.0.6" >> $job_script
  echo "module load libxml2/2.9.10" >> $job_script
  echo "module load R/4.2.1" >> $job_script
  echo "module load libx11-1.7.0-gcc-10.2.0-4iqvpeu" >> $job_script
  echo "module load libjpeg-turbo/2.0.6" >> $job_script
  echo "" >> $job_script
  echo "cd /scr/u/guhaogao/global_genomic_surveillance_spatpomp" >> $job_script
  echo "echo "Starting simulating ${MODEL} for the ${i} job."" >> $job_script
  echo "date" >> $job_script
  echo "Rscript scripts/model_simulation/Simulate_${MODEL}_2.R "results/model_data/model_simulation/${MODEL}/df_all_values_${MODEL}_${i}.rds" ${NCORES}" >> $job_script
  echo "date" >> $job_script
  echo "echo "Finished simulating ${MODEL} for the ${i} job"" >> $job_script
  echo "" >> $job_script
  echo "exit 0" >> $job_script

  sbatch $job_script
  rm $job_script
done
