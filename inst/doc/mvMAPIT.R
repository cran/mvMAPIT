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
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(kableExtra)
library(RcppAlgos)

## ----run_mapit_normal, eval = FALSE-------------------------------------------
#  mvmapit_hybrid <- mvmapit(
#          t(data$genotype),
#          t(data$trait),
#          test = "hybrid"
#  )
#  fisher <- fishers_combined(mvmapit_hybrid$pvalues)

## ----assign_data, include = FALSE---------------------------------------------
fisher <- mvmapit_data$fisher
pairs <- mvmapit_data$exhaustive_search
thresh <- 0.05 / (nrow(fisher) / 2)

## ----manhattan_plots----------------------------------------------------------
manhplot <- ggplot(fisher, aes(x = 1:nrow(fisher), y = -log10(p))) +
  geom_hline(yintercept = -log10(thresh), color = "grey40", linetype = "dashed") +
  geom_point(alpha = 0.75, color = "grey50") +
  geom_text_repel(aes(label=ifelse(p < thresh, as.character(id), '')))
plot(manhplot)

## ----significant_snps---------------------------------------------------------
thresh <- 0.05 / (nrow(fisher) / 2)
significant_snps <-  fisher %>%
    filter(p < thresh) # Call only marginally significant SNPs

truth <- simulated_data$epistatic %>%
  ungroup() %>%
  mutate(discovered = (name %in% significant_snps$id)) %>%
  select(name, discovered) %>%
  unique()

significant_snps %>%
  mutate_if(is.numeric, ~(as.character(signif(., 3)))) %>%
  mutate(true_pos = id %in% truth$name) %>%
  kable(., linesep = "") %>%
  kable_material(c("striped"))

## ----simulated_snps-----------------------------------------------------------
truth %>%
  kable(., linesep = "") %>%
  kable_material(c("striped"))

## ----search_significant_SNPs, eval = FALSE------------------------------------
#  # exhaustive search for p-values
#  pairs <- NULL
#  if (nrow(significant_snps) > 1) {
#    pairnames <- comboGeneral(significant_snps$id, 2)
#    # Generate unique pairs of SNP names;
#    # for length(names) = n, the result is a (n * (n-1)) x 2 matrix with one row corresponding to a pair
#    for (k in seq_len(nrow(pairnames))) {
#      fit <- lm(y ~ X[, pairnames[k, 1]]:X[, pairnames[k, 2]])
#      p_value1 <- coefficients(summary(fit))[[1]][2, 4]
#      p_value2 <- coefficients(summary(fit))[[2]][2, 4]
#      tib <- dplyr::tibble(
#              x = p_value1,
#              y = p_value2,
#              u = pairnames[k, 1],
#              v = pairnames[k, 2]
#      )
#      pairs <- bind_rows(pairs, tib)
#    }
#  }
#  
#  colnames(pairs) <- c(colnames(y), "var1", "var2")

## ----trait1_plots_exhaustive--------------------------------------------------
plotable <- pairs %>%
  pivot_longer(
    cols = starts_with("p_"),
    names_to = "trait",
    names_prefix = "trait ",
    values_to = "p",
    values_drop_na = TRUE
  ) %>%
  mutate(trait = case_when(
    trait == "p_01" ~ "Trait 1",
    trait == "p_02"  ~ "Trait 2"))
tiles <- ggplot(data = plotable, aes(x=var1, y=var2, fill=-log10(p)))+
  geom_tile() +
  facet_wrap(~trait) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_viridis_c()
plot(tiles)

## ----significant_exhaustive---------------------------------------------------
pairs %>%
  filter(p_01 < 0.05/15 | p_02 < 0.05/15) %>%
  kable(., linesep = "", digits = 14) %>%
  kable_material(c("striped"))

## ----simulated_interactions---------------------------------------------------
true_interactions <- simulated_data$interactions %>%
  mutate(var1 = sprintf(group1, fmt = "snp_%05d")) %>%
  mutate(var2 = sprintf(group2, fmt = "snp_%05d")) %>%
  mutate(trait = case_when(
    trait == 1 ~ "Trait 1",
    trait == 2  ~ "Trait 2")) %>%
  select(-c("group1", "group2"))
X <- true_interactions[, c("var1", "var2")]
X  <- t(apply(X, 1, sort))
true_interactions[,c("var1", "var2")]  <- X

epistatic_pairnames <- comboGeneral(simulated_data$epistatic$name %>% unique(), 2)
true_pairs <- NULL
for (k in seq_len(nrow(epistatic_pairnames))) {
  tib <- dplyr::tibble(var1 = epistatic_pairnames[k, 1],
                        var2 = epistatic_pairnames[k, 2])
  true_pairs <- bind_rows(true_pairs, tib)
}
anti <- anti_join(true_pairs, true_interactions) %>%
  mutate(effect_size = 0)

true_int_plot <- true_interactions %>%
  bind_rows(anti %>% mutate(trait = "Trait 1")) %>%
  bind_rows(anti %>% mutate(trait = "Trait 2"))

true_tiles <- ggplot(data = true_int_plot, aes(x=var1, y=var2, fill=effect_size)) +
  geom_tile() +
  facet_wrap(~trait) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white")
plot(true_tiles)

## ----true_interactions--------------------------------------------------------
true_interactions %>%
  kable(., linesep = "") %>%
  kable_material(c("striped"))

