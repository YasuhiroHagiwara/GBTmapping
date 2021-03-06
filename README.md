# GBT-Based Mapping Algorithms for EORTC QLQ-C30 onto EQ-5D-5L Index


1. This is a set of the R codes and RData files for mapping the EORTC QLQ-C30 onto the EQ-5D-5L index using regression-based mapping algorithms develped in "Hagiwara Y et al. Mapping EORTC QLQ-C30 and FACT-G onto EQ-5D-5L index for patients with cancer.Health Qual Life Outcomes. 20203;18:354." and gradient boosted tree (GBT)-based mapping algorithms developed in "Hagiwara Y et al. Gradient Boosted Tree Approaches for Mapping EORTC QLQ-C30 onto EQ-5D-5L Index for Patients with Cancer. Submitted."
2. If you use RStudio, you can use the R project file "R_functions_mapping.Rproj".
3. The main code is "R_functions_mapping.r", which includes the R functions for mapping.
4. In this R code, you can implement 4 mapping algorithms: direct and indirect mapping algorithms based on gradient boosted trees, direct mapping algorithms based on two-part beta regression, and indirect algorithms based on ordinal logistic regression.
5. Sample data are also provided, and sample EORTC QLQ-C30 data are mapped onto the EQ-5D-5L index in the provided R code.
