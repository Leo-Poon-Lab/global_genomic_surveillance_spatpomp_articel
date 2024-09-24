#!/bin/bash

if [ -z "$1" ]
  then
    echo "No model name supplied, use default model name: Omicron20"
    MODEL_NAME="Omicron20"
  else
    MODEL_NAME=$1
    # MODEL_NAME can only be "Omicron20" or "M1"
    if [ ${MODEL_NAME} != "Omicron20" ] && [ ${MODEL_NAME} != "M1" ]
      then
        echo "MODEL_NAME not valid!"
        exit 1
    fi
    echo "MODEL_NAME: ${MODEL_NAME}"
fi

## RUN ON LOCAL SERVER

echo "Start updating model"
date

if [ ${MODEL_NAME} == "Omicron20" ]
  then
    Rscript scripts/model_fitting/build_Omicron20.R
fi

if [ ${MODEL_NAME} == "M1" ]
  then
    Rscript scripts/model_fitting/build_M1.R
fi
date

echo "Finish updating model"
exit 0
