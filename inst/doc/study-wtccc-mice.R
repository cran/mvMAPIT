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

## ----run_mvmapit, eval = FALSE------------------------------------------------
#  mvmapit_TRAIT <- mvmapit(
#    t(TRAIT$genotype),
#    t(TRAIT$phenotype),
#    test = "hybrid"
#  )

## ----load_data, eval = FALSE--------------------------------------------------
#  mice_SI_paper <- readRDS("mice_SI_paper.rds")
#  mice_HCTHGB_MCVMCH <- readRDS("mice_HCTHGB_MCVMCH.rds")

## ----all_traits, eval = FALSE-------------------------------------------------
#  for_ticks_chr <- aggregate(position ~ chr, mice_data$fisher, function(x) c(first = min(x), last = max(x))) %>%
#    mutate(tick = floor((position[,"first"] + position[,"last"]) / 2)) %>%
#    mutate(chr2 = case_when(chr %% 5 == 0 ~ as.character(chr),
#                            chr == 1 ~ as.character(chr),
#                            TRUE ~ ""))
#  for_facetgrid_row <- as_labeller(c(`1` = "Trait #1", `2` = "Trait #2", `3` = "Covariance", `4` = "Combined"))
#  gg <- mice_SI_paper$fisher %>% ggplot(aes(
#        x = position,
#        y = -log10(pplot),
#        color = factor(color)
#      )) +
#        geom_point_rast(
#          size = 0.7) +
#        scale_color_manual(
#          values = c("#8b8b8b", "#bfbfbf", "#1b9e77")
#        ) +
#        scale_y_continuous(breaks = c(0, 5, 10),
#                           labels = c("0", "5", ">10")) +
#        geom_hline(
#          aes(
#            yintercept = -log10(5.179737e-06),
#            linetype = "Bonferroni"
#          ),
#          color = "#d95f02",
#          size = 0.3
#        ) +
#        theme_bw() +
#        facet_grid(x ~ y) +
#        theme(
#          panel.grid.major.x = element_blank(),
#          legend.position = "bottom",
#          text = element_text(family = "Times"),
#        ) +
#        labs(
#          y = bquote(-log[10](p)),
#          color = "") +
#        scale_x_continuous("Chromosome",
#                           breaks = for_ticks_chr$tick,
#                           labels = for_ticks_chr$chr2) +
#        scale_linetype_manual(name = "", values = c('dashed'))
#  show(gg)

## ----mice_data, eval = FALSE--------------------------------------------------
#  gg <- mice_HCTHGB_MCVMCH$fisher %>% ggplot(aes(
#        x = position,
#        y = -log10(pplot),
#        color = factor(color)
#      )) +
#        geom_point_rast(
#          size = 0.7) +
#        scale_color_manual(
#          values = c("#8b8b8b", "#bfbfbf", "#1b9e77")
#        ) +
#        scale_y_continuous(breaks = c(0, 5, 10),
#                           labels = c("0", "5", ">10")) +
#        geom_hline(
#          aes(
#            yintercept = -log10(5.179737e-06),
#            linetype = "Bonferroni"
#          ),
#          color = "#d95f02",
#          size = 0.3
#        ) +
#        theme_bw() +
#        facet_grid(row ~ case, labeller = labeller(row = for_facetgrid_row)) +
#        theme(
#          panel.grid.major.x = element_blank(),
#          legend.position = "bottom",
#          text = element_text(family = "Times"),
#        ) +
#        labs(
#          y = bquote(-log[10](p)),
#          color = "") +
#        scale_x_continuous("Chromosome",
#                           breaks = for_ticks_chr$tick,
#                           labels = for_ticks_chr$chr2) +
#        scale_linetype_manual(name = "", values = c('dashed'))
#  show(gg)

