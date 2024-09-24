#!/bin/bash

## RUN ON LOCAL SERVER

if [ -z "$1" ]
  then
    echo "INPUT_PATH not supplied!"
    exit 1
  else
    INPUT_PATH=$1
    echo "INPUT_PATH: ${INPUT_PATH}"
fi

if [ -z "$2" ]
  then
    echo "No model name supplied, use default model name: Omicron20"
    MODEL_NAME="Omicron20"
  else
    MODEL_NAME=$2
    # MODEL_NAME can only be "Omicron20" or "M1"
    if [ ${MODEL_NAME} != "Omicron20" ] && [ ${MODEL_NAME} != "M1" ]
      then
        echo "MODEL_NAME not valid!"
        exit 1
    fi
    echo "MODEL_NAME: ${MODEL_NAME}"
fi

echo "Starting mle extraction for ${INPUT_PATH}."
date
Rscript scripts/model_fitting/profiling_HKU/helper/analysis_local.R ${INPUT_PATH} ${MODEL_NAME}
date
echo "Finished mle extraction."

echo "Updating Figure 1."
date
Rscript scripts/model_simulation/Simulate_M0.R ${INPUT_PATH}
date
echo "Finished updating Figure 1."

exit 0