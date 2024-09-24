# Profiling and Local Searching on HKU HPC

Detailed commands and records can be found in [here](/draft/Model_fitting_steps.xlsx).

## Profiling and Local Searching for unit-specific/shared/all parameters

0. `0.update_model.sh` (on local server)
1. `1.profile_n.sh INPUT_MLE PARAM_TYPE ROUND_N PARTITION NCORES SEARCH_LEVEL ESTIMATOR MODEL_NAME` (on HPC)
2. `2.plot_profile_n.sh PARAM_TYPE ROUND_N MODEL_NAME` (on local server)
3. `3.local_n.sh INPUT_MLE LOCAL_PARAM_TYPE LOCAL_ROUND_N PARTITION NCORES SEARCH_LEVEL OPTIMIZER MODEL_NAME` (on HPC)
4. `4.plot_local_beta_n.sh INPUT_PATH MODEL_NAME` (on local server)

- Run step 0 to update the model frist.
- Run step 1 with `INPUT_MLE = NA` for the first time.
- `PARAM_TYPE` and `LOCAL_PARAM_TYPE` can be `unit`, `shared`, or `all`.
- Repeat steps 1, 2, 3, 4 in order with `ROUND_N` $\in 1...\infty$ .
- `PARTITION` can be "amd", "condo_amd" or "intel".
- `SEARCH_LEVEL` must be in "low", "mid", "high".
- `ESTIMATOR` can be "ubf" or "bpf".
- `OPTIMIZER` can be "iubf" or "ibpf".
- `MODEL_NAME` is optional, default is "Omicron20", can also be "M1".
- To make full use of the resources, if `SEARCH_LEVEL` is "low", `NCORES` over 32 is not allowed, if "mid", `NCORES` over 64 is not allowed, if "high", `NCORES` over 128 is not allowed.
- `LOCAL_ROUND_N` is the number of local searching rounds, can be run iteratively based on the previous local search result.
