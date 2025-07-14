# Dynamic Panel Data Selection Models

This project replicates the Monte Carlo Simulations from Section 3 of Sadoon et al. (2019).

**Note: This respository will be archived on 2025-08-14 or afterwards.**

## Instructions

### Replicating Table 1
1. Edit `code/run_parallel_N500.sh` and `code/run_parallel_N5000.sh` to use proper Stata path and navigate to `code/`.
2. Run `bash run_parallel_N500.sh`.
3. Run `bash run_parallel_N5000.sh`.
4. Run `python make_table_1.py`. You may have to activate the virtual environment or use `python3` instead of `python`.

### Replicating Table 2
1. Edit `code/run_experiment1.sh` through `code/run_experiment5.sh` to use proper Stata path and navigate to `code/`.
2. Run `bash run_experiment1.sh 500` through `bash run_experiment5.sh 500`.
3. Run `python make_table_2.py`. You may have to activate the virtual environment or use `python3` instead of `python`.

### Replicating Table 3
1. Edit `code/run_experiment1.sh` through `code/run_experiment5.sh` to use proper Stata path and navigate to `code/`.
2. Run `bash run_experiment1.sh 5000` through `bash run_experiment5.sh 5000`.
3. Run `python make_table_3.py`. You may have to activate the virtual environment or use `python3` instead of `python`.

### Replicating Figure 1
1. Edit `code/run_figure1_sims.sh` to use proper Stata path and navigate to `code/`.
2. Run `bash run_figure1_sims.sh` for all parameter combinations. The amount of simultaneous simulations will depend on the capabilities of your machine. For the small `N` values, I ran, for example:
   ```bash
   bash run_figure1_sims.sh "200 400 600 800 1000 1500 2000" "0.25" "A"
   ```
   But for larger `N` values I reduced the number of simultaneous simulations:
   ```bash
   bash run_figure1_sims.sh "2000 2500 3000" "0.25" "A"
   ```
3. Once all combinations have been run, execute `R fig1.R` to generate the figure. Output is three `.tex` files corresponding to TikZ figures.

