# GFLassoInference-experiment
Code and instructions to reproduce results and figures from the paper: 

Chen YT, Jewell SW, Witten DM. (2021) [More powerful selective inference for the graph fused lasso](https://arxiv.org/abs/2109.10451). arXiv:2109.10451 [statME].


### Step-by-step instructions for reproducing figures
#### Note that all the R codes assume that the working directory is `~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/`. Please double check your working directories and modify the relevant lines in the code accordingly (e.g., if you are running them on a computing cluster).

1. Download this repository and install relevant packages.
```
git clone https://github.com/yiqunchen/GFLassoInference-experiments.git
cd GFLassoInference-experiments
Rscript install_packages.R
```
2. Generate Figure 1 (Toy example to illustrate the higher power of our proposal).
```
Rscript ./Figure_1.R
```
3. Generate Figure 2 (Intuition for the perturbation).
```
Rscript ./Figure_2.R
```
4. Generate datasets for Figure 3 (1D fused lasso example).
```
Rscript ./Figure_3b_Gen_Data.R 
Rscript ./Figure_3c_Gen_Data.R 
```

5. Generate Figure 3, Figures 10(a)-(b), and Figure 11.
```
Rscript ./Figure_3a.R 
Rscript ./Figure_3b.R 
Rscript ./Figure_3c_11_10a_10b.R 
```
6. Generate datasets for Figure 4 (2D fused lasso example).
```
Rscript ./Figure_4b_Gen_Data.R
Rscript ./Figure_4c_Gen_Data.R
```
7. Generate Figure 4 (note that Figure 4(a) is the same as Figure 1(a)).

```
Rscript ./Figure_4b.R
Rscript ./Figure_4c_10c_10d.R
```

8. Generate Figure 5.

```
Rscript ./Figure_5.R
```

9. Generate Figure 6.

```
Rscript ./Figure_6.R
```

10. Generate Figure 8.
```
Rscript ./Figure_8a_Gen_Data.R
Rscript ./Figure_8b_Gen_Data.R
Rscript ./Figure_8.R
```
11. Generate Figure 9
```
Rscript ./Figure_9.R
```

12. Generate datasets for Figure 12
```
Rscript ./Figure_12_a_Gen_Data.R
Rscript ./Figure_12_b_Gen_Data.R
```
13. Generate Figure 12
```
Rscript ./Figure_12_a_b.R
```
14. Generate datasets for Figure 13
```
Rscript ./Figure_13_a_Gen_Data.R
Rscript ./Figure_13_b_Gen_Data.R
Rscript ./Figure_13_c_Gen_Data.R
```
15. Generate Figure 13
```
Rscript ./Figure_13_a.R
Rscript ./Figure_13_b.R
Rscript ./Figure_13_c.R
```
16. Generate Figure 14
```
Rscript ./Figure_14_Gen_Data.R
Rscript ./Figure_14.R
```
17. Generate Figure 15
```
Rscript ./Figure_15.R
```

18. Generate Figure 16
```
Rscript ./Figure_16.R
```

19. Generate Figure 17
```
Rscript ./Figure_17.R
```




