
####---- Specht et al. 2018 ---####


# Specht, Harrison, Guillaume Harmange, David H. Perlman, Edward Emmott, Zachary 
# Niziolek, Bogdan Budnik, and Nikolai Slavov. 2018. “Automated Sample 
# Preparation for High-Throughput Single-Cell Proteomics.” bioRxiv. 
# https://doi.org/10.1101/399774.

library(MSnbase)

# Download data
dat <- read.table("https://drive.google.com/uc?export=download&id=19o-vtyKOmVmlOoiPH4kri3mB4N3zUUcC", 
                  header = TRUE, sep = "\t")


