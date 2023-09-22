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
library(GGally)
library(tidyr)
library(dplyr)

## ----generate_data, eval = TRUE-----------------------------------------------
n_variants <- 10000
n_combine <- 3
pvalues <- tidyr::tibble(
   id = rep(as.character(c(1:n_variants)), each = n_combine),
   trait = rep(as.character(c(1:n_combine)), n_variants),
   p =  runif(n_variants * n_combine, min = 0, max = 1)
)

## ----combine, eval = TRUE-----------------------------------------------------
cauchy <- cauchy_combined(pvalues) %>%
  rename(p_cauchy = p) %>%
  select(-trait)
fisher <- fishers_combined(pvalues) %>%
  rename(p_fisher = p) %>%
  select(-trait)
harmonic <- harmonic_combined(pvalues) %>%
  rename(p_harmonic = p) %>%
  select(-trait)
min_max <- pvalues %>%
  group_by(id) %>%
  summarise(p_min = min(p), p_max = max(p))

combined_wide <- fisher %>%
  left_join(harmonic) %>%
  left_join(cauchy) %>%
  left_join(min_max) %>%
  select(-id)

## ----plot, eval = TRUE--------------------------------------------------------
my_bin <- function(data, mapping) {
  ggplot(data = data, mapping = mapping) +
  geom_bin2d() +
  scale_fill_continuous(type = "viridis")
}
ggpairs(combined_wide, columns = 1:5,
        lower = list(continuous = my_bin)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))

