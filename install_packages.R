packages <- c("devtools", "tidyverse", "latex2exp", "RColorBrewer",
              "shades", "ggpubr","igraph", "genlasso")
install.packages(setdiff(packages, rownames(installed.packages())))  
if(!("GFLassoInference" %in% rownames(installed.packages()))){
  devtools::install_github("yiqunchen/GFLassoInference")
}
