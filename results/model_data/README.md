## Key model files

| File path                          | Description                                                                                                                              |
|---------------------------|---------------------------------------------|
| profiling3_Omicron20/mlev2_mid.rds | The highest log-likelihood (loglik) parameter set we have achieved so far, representing the current Maximum Likelihood Estimation (MLE). |
| profile3_Omicron20_local1          | A local search based on profiling3_Omicron20/mle_mid.rds, which did not further improve the loglik.                                      |
| profile3_Omicron20_local1v2        | Another local search based on profiling3_Omicron20/mlev2_mid.rds, also without further loglik improvement.                               |

## Sequence of model developments

Below is a representation of the relationship and sequence of the model files:

``` mermaid
graph TD
    A[profiling1_Omicron20/mle_mid.rds]
    B[profiling1_Omicron20_local1/mle_mid.rds]
    C[profiling2_Omicron20/mlev2_mid.rds]
    D[profiling2_Omicron20_local1v2/mle_mid.rds]
    E[profiling3_Omicron20/mlev2_mid.rds]

    A --> B
    B --> C
    C --> D
    D --> E
```
