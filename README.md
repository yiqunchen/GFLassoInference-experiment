# GFLassoInference-experiment
Code and instructions to reproduce results and figures from the paper: 

Chen YT, Jewell SW, Witten DM. (2021) [More powerful selective inference for the graph fused lasso](). arXiv:xxxx.xxxx [statME].


### Step-by-step instructions for reproducing figures
#### Note that all the R codes assume that the working directory is `~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/`. Please double check your working directories and modify the relevant lines in the code accordingly (e.g., if you are running them on a computing cluster).

1. Download this repository and install relevant packages.
```
git clone https://github.com/yiqunchen/GFLassoInference-experiments.git
cd GFLassoInference-experiments
Rscript install_packages.R
```
2. Generate Figure 1 (Toy example for illustrating our proposal).
```
Rscript ./Figure_1.R
```
3. Generate Figure 2 (Intuition for perturbation dataset).
```
Rscript ./Figure_2.R
```
4. Generate datasets for Figure 3 (1D fused lasso example).
```
Rscript ./Figure_3b_Gen_Data.R 
Rscript ./Figure_3c_Gen_Data.R 
```

5. Generate Figure 3 and Figure 7(a).
```
Rscript ./Figure_3a.R 
Rscript ./Figure_3b.R 
Rscript ./Figure_3c_7a.R 
```
6. Generate datasets for Figure 4 (2D fused lasso example).
```
Rscript ./Figure_4b_Gen_Data.R
Rscript ./Figure_4c_Gen_Data.R
```
7. Generate Figure 4 (note that Figure 4(a) is the same as Figure 1(a)).
```
Rscript ./Figure_4b.R
Rscript ./Figure_4c.R
```
8. Generate Figure 5.
```
Rscript ./Figure_5.R
```
9. Generate Figure 6.
```
Rscript ./Figure_6.R
```
11. Generate Figure 7 (note that Figure 7(a) is generated in step 5).

```
Rscript ./Figure_7b_Gen_Data.R
Rscript ./Figure_7b.R
```

12. Generate datasets for Figure 8

```
Rscript ./Figure_8a_Gen_Data.R
Rscript ./Figure_8b_Gen_Data.R
```

13. Generate Figure 8

```
Rscript ./Figure_8a_b.R
```





