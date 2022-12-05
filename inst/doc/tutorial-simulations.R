## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#",
  fig.width=7,
  fig.height=5
)
library(knitr)

## ----load_dependencies, message=FALSE-----------------------------------------
library(mvMAPIT)

## ----simulate_genotypes, eval = FALSE-----------------------------------------
#  set.seed(1234)
#  
#  n_samples <- 2938
#  n_snp <- 5747
#  
#  maf <- 0.05 + 0.45 * runif(n_snp)
#  random_minor_allele_counts   <- (runif(n_samples * n_snp) < maf) + (runif(n_samples * n_snp) < maf)
#  genotype_data <- matrix(random_minor_allele_counts,
#    nrow = n_samples,
#    ncol = n_snp,
#    byrow = TRUE,
#  )
#  
#  sample_names <- seq_len(n_samples) %>% sprintf(fmt = "id%04d")
#  snp_names <- seq_len(n_snp) %>% sprintf(fmt = "snp%04d")
#  
#  colnames(genotype_data) <- snp_names
#  rownames(genotype_data) <- sample_names

## ----simulation_parameters, eval = FALSE--------------------------------------
#  seed <- 67132
#  d <- 2
#  PVE <- 0.6
#  rho <- 0.2
#  n_causal <- 1000
#  n_trait_specific <- 0
#  n_pleiotropic <- 10
#  group_ratio_pleiotropic <- 1
#  epistatic_correlation <- 0.9
#  maf <- 0.05
#  simulated_data <- simulate_traits(
#      genotype_data,
#      n_causal = n_causal,
#      n_trait_specific = n_trait_specific,
#      n_pleiotropic = n_pleiotropic,
#      d = d,
#      H2 = PVE,
#      rho = rho,
#      epistatic_correlation = epistatic_correlation,
#      group_ratio_pleiotropic = group_ratio_pleiotropic,
#      maf_threshold = maf,
#      seed = seed
#  )

