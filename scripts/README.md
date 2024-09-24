The workflow follows the below order:
1. [data_preprocessing](data_processing/)
    - [`install_prerequisites.sh`](data_processing/install_prerequisites.sh)
    - [`Prepare_IHME_data.r`](data_processing/Prepare_IHME_data.r)
    The IHME mainly includes the 1. daily reported cases (raw data), and 2. daily reported deaths (raw data with excess mortality scalar applied), for various spatial units and is used for the model fitting. We also used the below data from IHME as covariates while model fitting: 1. Infection detection ratio, 2. Infection fatality ratio, 3. Effectively vaccinated individuals (one and two dose with efficacy). 
    - [`Prepare_GISAID_data.r`](data_processing/Prepare_GISAID_data.r)
    The GISAID data represents the daily sequenced cases for various viral lineages (as Pango lineages, which will later be assigned to our studied strains Delta, Omicron BA.1 and Omicron BA.2 variants) and is used for the model fitting.
    - [`Prepare_origin_data.r`](data_processing/Prepare_origin_data.r)
    We infer the origin of different sequences using phylogenetic methods. This data is used for the model fitting.
    - [`Prepare_travel_data.r`](data_processing/Prepare_travel_data.r)
    This script is for analysing and generating the global movement matrix between regions, basing on three different datasets. Details please refer to the methods section of the manuscript.
    - [`infer_origins`](data_processing/infer_origins/)
    This folder contains the scripts for inferring the origin of the sequences. The script [`infer_origins.sh`](data_processing/infer_origins/infer_origins.sh) is used to run the inference on the cluster.

2. [model_building](model_building/)
    - The script [`spatPomp_Omicron20.R`](model_building/spatPomp_Omicron20.R) contains a function used for building model $M_0$. Similarly, [`spatPomp_M1.R`](model_building/spatPomp_M1.R), [`spatPomp_M3.R`](model_building/spatPomp_M3.R), [`spatPomp_M4.R`](model_building/spatPomp_M4.R) contain functions for building models $M_1$, $M_3$, and $M_4$, respectively. 

3. [model_fitting](model_fitting/)
    - Modeling fitting scripts for the *burn-in* and the formal fitting of models $M_0$ and $M_1$ are available in the folder [`profiling_HKU`](model_fitting/profiling_HKU/). 
    - Specification and fitting of the model $M_2$ are available in the script [`benchmark.R`](model_fitting/benchmark.R).

    Details are provided in the subfolder.

4. [model_simulation](model_simulation/)
    This folder contains the scripts for simulating the models $M_0$, $M_3$, and $M_4$. Details are provided in the subfolder.

