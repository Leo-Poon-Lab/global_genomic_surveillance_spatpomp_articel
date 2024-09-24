## Simulating $M_0$

- Run `Rscript scripts/model_simulation/Simulate_M0.R MLE_PATH` with the path of input_mle.

## Simulating $M_3$

1. Run `Rscript scripts/model_simulation/Simulate_M3_1.R` to update the `Simulate_M3_2_HPC_command.sh` file.
2. **Manually update** the `path_mle` in `scripts/model_simulation/Simulate_M3_2.R` to the path of the latest MLE.
3. Run `scripts/model_simulation/Simulate_M3_2_HPC_command.sh` **on HPC** to simulate $M_3$.
4. Run `Rscript scripts/model_simulation/Simulate_M3_3.R` to update the plot.

## Simulating $M_4$

1. Run `Rscript scripts/model_simulation/Simulate_M4_1.R` to update the `Simulate_M4_2_HPC_command.sh` file.
2. **Manually update** the `path_mle` in `scripts/model_simulation/Simulate_M4_2.R` to the path of the latest MLE.
3. Run `scripts/model_simulation/Simulate_M4_2_HPC_command.sh` **on HPC** to simulate $M_4$.
4. Run `Rscript scripts/model_simulation/Simulate_M4_3.R` to update the plot and generate data for the second simulation round.
5. Run `scripts/model_simulation/Simulate_M4_3_HPC_command.sh` **on HPC** to simulate $M_4$ for the second round.
6. Run `Rscript scripts/model_simulation/Simulate_M4_4.R` to update the plot.
