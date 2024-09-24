#!/bin/bash

if [ -z "$1" ]
  then
    echo "INPUT_MLE not supplied!"
    exit 1
  else
    INPUT_MLE=$1
    echo "INPUT_MLE: ${INPUT_MLE}"
fi

if [ -z "$2" ]
  then
    echo "PARAM_TYPE not supplied!"
    exit 1
  else
    PARAM_TYPE=$2
    echo "PARAM_TYPE: ${PARAM_TYPE}"
fi
# PARAM_TYPE must be in "all" or "unit" or "shared"
if [ ${PARAM_TYPE} != "all" ] && [ ${PARAM_TYPE} != "unit" ] && [ ${PARAM_TYPE} != "shared" ]
  then
    echo "PARAM_TYPE not valid!"
    exit 1
fi
# If PARAM_TYPE is "all", profile all parameters, set N_PARAMS to 49, if "unit", profile unit parameters, set N_PARAMS to 29, if "shared", profile shared parameters, set N_PARAMS to 20
if [ ${PARAM_TYPE} == "all" ]
  then
    N_PARAMS=49
elif [ ${PARAM_TYPE} == "unit" ]
  then
    N_PARAMS=29
elif [ ${PARAM_TYPE} == "shared" ]
  then
    N_PARAMS=20
fi

if [ -z "$3" ]
  then
    echo "ROUND_N not supplied!"
    exit 1
  else
    ROUND_N=$3
    echo "ROUND_N: ${ROUND_N}"
fi

if [ -z "$4" ]
  then
    echo "PARTITION not supplied!"
    exit 1
  else
    PARTITION=$4
    echo "PARTITION: ${PARTITION}"
fi
# PARTITIONS must be in "amd", "condo_amd" or "intel"
if [ ${PARTITION} != "amd" ] && [ ${PARTITION} != "condo_amd" ] && [ ${PARTITION} != "intel" ]
  then
    echo "PARTITION not valid!"
    exit 1
fi
# profile every parameters
if [ ${PARTITION} == "amd" ]
  then
  seq_i=$(seq ${N_PARAMS} -1 1)
  # seq_i=$(seq 17 -1 8)
else
  seq_i=$(seq 1 ${N_PARAMS})
  # seq_i=$(seq 10 11)
fi

if [ -z "$5" ]
  then
    echo "NCORES not supplied!"
    exit 1
  else
    NCORES=$5
    echo "NCORES: ${NCORES}"
fi

if [ -z "$6" ]
  then
    echo "No SEARCH_LEVEL supplied"
    exit 1
  else
    SEARCH_LEVEL=$6
    echo "SEARCH_LEVEL: ${SEARCH_LEVEL}"
fi
# SEARCH_LEVEL must be in "low", "mid", "high"
if [ ${SEARCH_LEVEL} != "low" ] && [ ${SEARCH_LEVEL} != "mid" ] && [ ${SEARCH_LEVEL} != "high" ]
  then
    echo "SEARCH_LEVEL not valid!"
    exit 1
fi
# if SEARCH_LEVEL is "low", NCORES over 32 is not allowed, if "mid", NCORES over 64 is not allowed, if "high", NCORES over 128 is not allowed
if [ ${SEARCH_LEVEL} == "low" ] && [ ${NCORES} -gt 32 ]
  then
    echo "For low search level, NCORES over 32 is not recommended!"
    exit 1
fi
if [ ${SEARCH_LEVEL} == "mid" ] && [ ${NCORES} -gt 64 ]
  then
    echo "For mid search level, NCORES over 64 is not recommended!"
    exit 1
fi
if [ ${SEARCH_LEVEL} == "high" ] && [ ${NCORES} -gt 128 ]
  then
    echo "For high search level, NCORES over 128 is not recommended!"
    exit 1
fi

if [ -z "$7" ]
  then
    echo "No estimator algorithm supplied, use ubf as defualt"
    ESTIMATOR="ubf"
  else
    ESTIMATOR=$7
    echo "ESTIMATOR: ${ESTIMATOR}"
fi
# ESTIMATOR must be either "ubf" or "bpf"
if [ ${ESTIMATOR} != "ubf" ] && [ ${ESTIMATOR} != "bpf" ]
  then
    echo "ESTIMATOR not valid!"
    exit 1
fi
# if ESTIMATOR is "ubf", memory per CPU is 1200M, if "bpf", memory is 3200M
if [ ${ESTIMATOR} == "ubf" ]
  then
    MEM_PER_CPU="1200M"
elif [ ${ESTIMATOR} == "bpf" ]
  then
    MEM_PER_CPU="3200M"
fi

if [ -z "$8" ]
  then
    echo "No model name supplied, use default model name: Omicron20"
    MODEL_NAME="Omicron20"
  else
    MODEL_NAME=$8
    # MODEL_NAME can only be "Omicron20" or "M1"
    if [ ${MODEL_NAME} != "Omicron20" ] && [ ${MODEL_NAME} != "M1" ]
      then
        echo "MODEL_NAME not valid!"
        exit 1
    fi
    echo "MODEL_NAME: ${MODEL_NAME}"
fi

if [ ${MODEL_NAME} == "M1" ]
  then
    if [ ${SEARCH_LEVEL} == "low" ]
      then
        TIME="00:30:00"
    elif [ ${SEARCH_LEVEL} == "mid" ]
      then
        TIME="02:00:00"
    elif [ ${SEARCH_LEVEL} == "high" ]
      then
        TIME="05:00:00"
    fi
  else
    if [ ${SEARCH_LEVEL} == "low" ]
      then
        TIME="7:00:00"
    elif [ ${SEARCH_LEVEL} == "mid" ]
      then
        TIME="15:00:00"
    elif [ ${SEARCH_LEVEL} == "high" ]
      then
        if [ ${NCORES} -lt 64 ]
        then
          TIME="48:00:00"
        else
          TIME="24:00:00"
        fi
    fi
fi

## RUN ON HPC
for i in $seq_i; do
  job_script="job_$i.sh"
  echo "#!/bin/bash" > $job_script
  echo "" >> $job_script
  echo "## NECESSARY JOB SPECIFICATIONS" >> $job_script
  echo "#SBATCH --job-name=p_${PARAM_TYPE}_${ROUND_N}_$i" >> $job_script
  echo "#SBATCH --partition=${PARTITION}" >> $job_script
  echo "#SBATCH --qos=normal" >> $job_script
  echo "#SBATCH --time=${TIME}" >> $job_script
  echo "#SBATCH --nodes=1" >> $job_script
  echo "#SBATCH --ntasks=1" >> $job_script
  echo "#SBATCH --cpus-per-task=${NCORES}" >> $job_script
  echo "#SBATCH --exclude=GPA-2-9" >> $job_script
  echo "#SBATCH --mem-per-cpu=${MEM_PER_CPU}" >> $job_script
  echo "#SBATCH --output=%x_%j.out" >> $job_script
  echo "" >> $job_script
  echo "## OPTIONAL JOB SPECIFICATIONS" >> $job_script
  echo "#SBATCH --mail-user=guhaogao@hku.hk" >> $job_script
  echo "#SBATCH --mail-type=END,FAIL" >> $job_script
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
  echo "date" >> $job_script
  echo "TASK_ID=$i Rscript scripts/model_fitting/profiling_HKU/helper/profile_params_n.R ${INPUT_MLE} ${PARAM_TYPE} ${ROUND_N} ${NCORES} ${SEARCH_LEVEL} ${ESTIMATOR} ${MODEL_NAME}" >> $job_script  
  echo "date" >> $job_script
  echo "exit 0" >> $job_script

  sbatch $job_script
  rm $job_script
done
