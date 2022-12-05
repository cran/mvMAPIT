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
library(ggplot2)
library(dplyr)

## ----read_data, eval = FALSE--------------------------------------------------
#  data.files <- c("elife-71393-fig1-data1-v3.csv", "elife-71393-fig1-data2-v3.csv")
#  
#  data1 <- read.csv(file = data.files[1])
#  data2 <- read.csv(file = data.files[2])
#  
#  data1names <- c(30, 35, 36, 57, 64, 65, 66, 79, 82, 83, 84, 85, 92, 95, 103, 113)
#  data2names <- c(29, 35, 65, 66, 69, 82, 83, 84, 85, 87, 112.1)
#  data1map <- as.data.frame(list(
#    pos = sprintf("pos%d", c(1:16)),
#    term = data1names
#  ))
#  data2map <- as.data.frame(list(
#    pos = sprintf("pos%d", c(1:11)),
#    term = data2names
#  ))
#  phenotype1 <- data1[, endsWith(names(data1), "mean")]
#  genotype1 <- data1[, startsWith(names(data1), "pos")]
#  
#  phenotype2 <- data2[, endsWith(names(data2), "mean")]
#  genotype2 <- data2[, startsWith(names(data2), "pos")]
#  
#  colnames(genotype1) <- data1names
#  colnames(genotype2) <- data2names
#  
#  CR9114 <- list("phenotype" = phenotype1, "genotype" = genotype1, "map" = data1map)
#  CR6261 <- list("phenotype" = phenotype2, "genotype" = genotype2, "map" = data2map)

## ----run_mvmapit, eval = FALSE------------------------------------------------
#  mvmapit_CR9114 <- mvmapit(
#    t(CR9114$genotype),
#    t(CR9114$phenotype),
#    test = "hybrid"
#  )
#  mvmapit_CR6261 <- mvmapit(
#    t(CR6261$genotype),
#    t(CR6261$phenotype),
#    test = "hybrid"
#  )

## ----manhattan----------------------------------------------------------------
for_facetgrid_row <-
  as_labeller(c(
    `1` = "Trait #1",
    `2` = "Trait #2",
    `3` = "Covariance",
    `4` = "Combined"
  ))
phillips_data$fisher$colorf <- factor(phillips_data$fisher$color,
                                      labels = c("1",
                                                 "Significant"))
gg_fisher <- phillips_data$fisher %>%
  ggplot(aes(x = position, y = -log10(pplot))) +
  geom_point(aes(colour = colorf), size = 1) +
  scale_color_manual(values = c("#1b9e77", "#2c2c2c"),
                     breaks = c("Significant")) +
  scale_y_continuous(breaks = c(0, 5, 10),
                     labels = c("0", "5", ">10")) +
  geom_hline(aes(yintercept = -log10(threshold),
                 linetype = "Bonferroni"),
             color = "#d95f02") +
  scale_linetype_manual(name = "", values = c('dashed')) +
  theme_bw() +
  facet_grid(row ~ species, labeller = labeller(row = for_facetgrid_row)) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    text = element_text(family = "Arial")
  ) +
  labs(x = "Position",
       y = "-log10(p)",
       colour = NULL)
show(gg_fisher)

## ----CR6261-------------------------------------------------------------------
gg_CR6261 <- ggplot(phillips_data$regression$CR6261, 
                    aes(res_x, res_y, fill = effect)) + 
  geom_tile() +
  scale_fill_gradient2(
    high = "#b2182b",
    mid = "white",
    low = "#2166ac",
    midpoint = 0,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")
show(gg_CR6261)

## ----CR9114-------------------------------------------------------------------
gg_CR9114 <- ggplot(phillips_data$regression$CR9114, 
                    aes(res_x, res_y, fill = effect)) + 
  geom_tile() +
  scale_fill_gradient2(
    high = "#b2182b",
    mid = "white",
    low = "#2166ac",
    midpoint = 0,
    space = "Lab",
    na.value = "#000000",
    guide = "colourbar",
    limits = c(min(phillips_data$regression$CR6261$effect), 
               max(phillips_data$regression$CR6261$effect))
  ) +
  theme_bw() +
  theme(legend.position = "bottom")
show(gg_CR9114)

