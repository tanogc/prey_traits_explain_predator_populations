# Calculating prey biodiversity metrics to predict propulations of high-value predators

Code to re-create analysis for "Populations of high-value predators reflect the traits of their prey". 

In this study, we characterise prey use and distribution in iconic bird (grey wagtails and Eurasian dippers) and fish species (brown trout and Atlantic salmon) to assess whether prey traits could predict populations of these four riverine predators. Specifically, we hypothesised that: (i) Prey key traits would predict predator populations more effectively than (ii) Diversity of prey traits, (iii) the taxonomic abundance or richness of prey (known as traditional or mass-effect types of biodiversity) or (iv) the prevailing environmental conditions. 

This code reproduces analysis to estimate diet preferences of the four riverine predators, model exploring the predicting capacity of the four hypothesised mechanisms and a multi-threshold analysis to determine to what extent these four predator species could be sustained simultaneously under current environmental conditions.

## Original article:

Please, use this citation to reference the code:

```
Gutiérrez-Cánovas, C., Worthington, T.A., Jâms, I.B., Noble, D.G., Perkins, D.M., 
Vaughan, I.P., Woodward, G., Ormerod, S.J. & I. Durance, 2020. Populations of 
high-value predators reflect the traits of their prey. Ecography.
```

## R files description:

* **0_FD_functions.R**: R script to estimate Functional Diversity (FD) metrics
* **0_quality_funct_space_fromdist.R**: R function for computing the quality of functional dendrogramm and multidimensional functional spaces. This function is a simplified version of the Appendix S1 associated to Maire et al. 2015 (Global Ecol. and Biogeogr.)
* **1_functional_groups.R**: R script to classify invertebrate prey into functional groups based on effect traits
* **2_diet analysis.R**: R script to analyse diet preferences of riverine predators
* **3_preparing_data.R**: R script to quantify prey biodiversity metrics and prepare response variables (predator populations)
* **4_models.R**: R script to relate biodiversity dimensions with single and multiple predator population sizes
* **5_multi_threshold.R**: R script to perform the multi-threshold analysis

## Original data
* **aggregated_diet.txt**: data on diet preferences extracted from literature of the four riverine predators
* **e_traits.txt**: effect traits for the taxa occuring in the study area (89 taxa)
* **e_traits_all.txt**: effect traits for the taxa occuring in the diet preference data from literature (175 taxa)
* **r_traits.txt**: effect traits for the taxa occuring in the study area (89 taxa)
* **r_traits_all.txt**: effect traits for the taxa occuring in the diet preference data from literature (175 taxa)

## Dependencies
To run the code and functions from this repository, you need to install the following packages: "ape", "arm", "car", "ecodist", "FD", "lme4", "lmerTest", 
"multcomp", "multifunc", "MuMIn", "nlme", "plyr", "reshape2", "sandwich", "sqldf", "sqldf", "usdm", "vegan", "viridis". Use this code to install them:

```
install.packages(c("ape", "arm", "car", "ecodist", "FD", "lme4", "lmerTest", 
"multcomp", "MuMIn", "nlme", "plyr", "reshape2", "sandwich", "sqldf", "sqldf", 
"usdm", "vegan", "viridis")
              
library(devtools)
install_github("jebyrnes/multifunc", force=T)
```

Please, send questions or problems related with the use of this code to Cayetano Gutiérrez Cánovas (cayeguti@um.es).

