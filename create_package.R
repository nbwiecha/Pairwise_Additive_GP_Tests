# Create R package to store all compiled C++ code

library(Rcpp)

setwd("~/GitHub/Pairwise_Additive_GP_Tests")

# Rcpp.package.skeleton("PairwiseSims")
devtools::clean_dll("PairwiseSims")
compileAttributes("PairwiseSims")
# install.packages("devtools") # If you don't already have it
devtools::build("PairwiseSims")
# devtools::document("PairwiseSims")
install.packages("PairwiseSims_0.0.1.tar.gz", repos=NULL, source=TRUE)
