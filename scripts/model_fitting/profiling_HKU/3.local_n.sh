#!/bin/bash

## RUN ON HPC

if [ -z "$1" ]
  then
    echo "INPUT_MLE not supplied!"
    exit 1
  else
    INPUT_MLE=$1
    INPUT_PATH=$(dirname ${INPUT_MLE})
    INPUT_SEARCH=$(echo ${INPUT_PATH} | sed 's/.*profiling_//g;')
    echo "INPUT_MLE: ${INPUT_MLE}"
fi

if [ -z "$2" ]
  then
    echo "LOCAL_PARAM_TYPE not supplied!"
    exit 1
  else
    LOCAL_PARAM_TYPE=$2
    echo "LOCAL_PARAM_TYPE: ${LOCAL_PARAM_TYPE}"
fi
# LOCAL_PARAM_TYPE must be in "all" or "unit" or "shared"
if [ ${LOCAL_PARAM_TYPE} != "all" ] && [ ${LOCAL_PARAM_TYPE} != "unit" ] && [ ${LOCAL_PARAM_TYPE} != "shared" ]
  then
    echo "LOCAL_PARAM_TYPE not valid!"
    exit 1
fi

if [ -z "$3" ]
  then
    echo "LOCAL_ROUND_N not supplied!"
    exit 1
  else
    LOCAL_ROUND_N=$3
    OUTPUT_PATH=$(echo ${INPUT_PATH} | sed -e 's#_local_[0-9]\+##g;')_local_${LOCAL_ROUND_N}
    echo "LOCAL_ROUND_N: ${LOCAL_ROUND_N}"
    echo "OUTPUT_PATH: ${OUTPUT_PATH}"
fi

if [ -z "$4" ]
  then
    echo "PARTITION not supplied!"
    exit 1
  else
    PARTITION=$4
    echo "PARTITION: ${PARTITION}"
fi
# PARTITIONS must be in "amd", "condo_amd" or "intel" or "hugemem"
if [ ${PARTITION} != "amd" ] && [ ${PARTITION} != "condo_amd" ] && [ ${PARTITION} != "intel" ] && [ ${PARTITION} != "hugemem" ]
  then
    echo "PARTITION not valid!"
    exit 1
fi
if [ ${PARTITION} == "hugemem" ]
  then
    QOS="hugemem"
  else
    QOS="normal"
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
    echo "No search_level supplied"
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
    echo "No optimizing algorithm supplied, use iubf as defualt"
    OPTIMIZER="iubf"
  else
    OPTIMIZER=$7
    echo "OPTIMIZER: ${OPTIMIZER}"
fi
# OPTIMIZER must be either "iubf" or "ibpf"
if [ ${OPTIMIZER} != "iubf" ] && [ ${OPTIMIZER} != "ibpf" ]
  then
    echo "OPTIMIZER not valid!"
    exit 1
fi
# if OPTIMIZER is "iubf", memory per CPU is 1200M, if "ibpf", memory is 3200M
if [ ${OPTIMIZER} == "iubf" ]
  then
    MEM_PER_CPU="1200M"
elif [ ${OPTIMIZER} == "ibpf" ]
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
    elif [ ${SEARCH_LEVEL} == "high" ] # nubf=1
      then
      if [ ${NCORES} -lt 64 ]
        then
          TIME="5:00:00"
        else
          TIME="3:00:00"
        fi
    fi
fi

# N_nodes
N_nodes=20

for i in $(seq 1 ${N_nodes}); do
  job_script="job_$i.sh"
  echo "#!/bin/bash" > $job_script
  echo "" >> $job_script
  echo "## NECESSARY JOB SPECIFICATIONS" >> $job_script
  echo "#SBATCH --job-name=$(basename ${OUTPUT_PATH} | sed 's/profiling_//g;')_${LOCAL_PARAM_TYPE}_${PARTITION}_${NCORES}_${i}" >> $job_script
  echo "#SBATCH --partition=${PARTITION}" >> $job_script
  echo "#SBATCH --qos=${QOS}" >> $job_script
  echo "#SBATCH --time=${TIME}" >> $job_script
  echo "#SBATCH --nodes=1" >> $job_script
  echo "#SBATCH --ntasks=1" >> $job_script
  echo "#SBATCH --cpus-per-task=${NCORES}" >> $job_script
  echo "#SBATCH --mem-per-cpu=${MEM_PER_CPU}" >> $job_script
  echo "#SBATCH --output=%x_%j.out" >> $job_script
  echo "" >> $job_script
  echo "## OPTIONAL JOB SPECIFICATIONS" >> $job_script
  echo "#SBATCH --mail-user=guhaogao@hku.hk" >> $job_script
  echo "#SBATCH --mail-type=END,FAIL" >> $job_script
  echo "" >> $job_script
  echo "" >> $job_script
  echo "module load gcc/10.2" >> $job_script
  echo "module load cmake/3.19.7" >> $job_script
  echo "module load libjpeg-turbo/2.0.6" >> $job_script
  echo "module load libxml2/2.9.10" >> $job_script
  echo "module load R/4.2.1" >> $job_script
  echo "module load libx11-1.7.0-gcc-10.2.0-4iqvpeu" >> $job_script
  echo "module load libjpeg-turbo/2.0.6" >> $job_script
  echo "" >> $job_script
  #echo "cd /scr/u/guhaogao/global_genomic_surveillance_spatpomp" >> $job_script
  echo "echo "Starting ${PARAM_TYPE} local search for ${INPUT_SEARCH} ${ROUND_N}_${i}."" >> $job_script
  echo "date" >> $job_script
  
  echo "    echo "using ${OPTIMIZER} for ${INPUT_MLE}."" >> $job_script
  echo "    Rscript scripts/model_fitting/profiling_HKU/helper/local_params_n.R ${INPUT_MLE} ${LOCAL_PARAM_TYPE} ${OUTPUT_PATH}/ ${SEARCH_LEVEL} ${PARTITION} ${NCORES} ${i} ${OPTIMIZER} ${MODEL_NAME}" >> $job_script

  echo "date" >> $job_script
  echo "echo "Finished local search"" >> $job_script
  echo "" >> $job_script
  echo "exit 0" >> $job_script

  sbatch $job_script
  rm $job_script
done
