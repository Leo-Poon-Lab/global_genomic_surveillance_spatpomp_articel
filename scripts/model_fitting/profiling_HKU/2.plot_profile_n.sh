#!/bin/bash

## RUN ON LOCAL SERVER

if [ -z "$1" ]
  then
    echo "PARAM_TYPE not supplied!"
    exit 1
  else
    PARAM_TYPE=$1
    echo "PARAM_TYPE: ${PARAM_TYPE}"
fi
# PARAM_TYPE must be in "all" or "unit" or "shared"
if [ ${PARAM_TYPE} != "all" ] && [ ${PARAM_TYPE} != "unit" ] && [ ${PARAM_TYPE} != "shared" ]
  then
    echo "PARAM_TYPE not valid!"
    exit 1
fi

if [ -z "$2" ]
  then
    echo "ROUND_N not supplied!"
    exit 1
  else
    ROUND_N=$2
    echo "ROUND_N: ${ROUND_N}"
fi

if [ -z "$3" ]
  then
    echo "No model name supplied, use default model name: Omicron20"
    MODEL_NAME="Omicron20"
  else
    MODEL_NAME=$3
    # MODEL_NAME can only be "Omicron20" or "M1"
    if [ ${MODEL_NAME} != "Omicron20" ] && [ ${MODEL_NAME} != "M1" ]
      then
        echo "MODEL_NAME not valid!"
        exit 1
    fi
    echo "MODEL_NAME: ${MODEL_NAME}"
fi

echo "Starting plotting."
date
Rscript scripts/model_fitting/profiling_HKU/helper/plot_profile.R "results/model_data/profiling_HKU/${MODEL_NAME}/profiling_${PARAM_TYPE}_${ROUND_N}" "results/figs/profiling_HKU/${MODEL_NAME}/profiling_${PARAM_TYPE}_${ROUND_N}"
date
echo "Finished plotting."

echo "Starting mle extraction."
date
Rscript scripts/model_fitting/profiling_HKU/helper/analysis_profile.R "results/model_data/profiling_HKU/${MODEL_NAME}/profiling_${PARAM_TYPE}_${ROUND_N}" ${MODEL_NAME} 
date
echo "Finished mle extraction."

echo "Updating Figure 1."
date
Rscript scripts/model_simulation/Simulate_M0.R "results/model_data/profiling_HKU/${MODEL_NAME}/profiling_${PARAM_TYPE}_${ROUND_N}/"
date
echo "Finished updating Figure 1."

exit 0