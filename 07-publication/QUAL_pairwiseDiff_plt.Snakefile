################################################################################
# Project: Ancient human remains from Guam
# Part: Preparing material for publication
# Step: Plot for pairwise difference analysis
#
# Alex Huebner, 26/03/2020
################################################################################

from snakemake.utils import R

workdir: "../snakemake_tmp"

rule all:
    input: 


rule pairwiseDiff_plt:
    output:
        "supp_figures/QUAL_pairwiseDist.pdf"
    message: "Plot the pairwise difference analysis results for the Guam samples"
    params: 
        dist = "../04-analysis/qual/pairwiseDist/pairwiseDist_1240K.txt.gz"
    run:
        R("""
          library(data.table)
          library(tidyverse)
          
          pairwiseDiff <- fread("{params.dist}")
          
          # Extract all comparisons with the Guam samples
          guam_comparisons <- pairwiseDiff %>%
                              filter(ind1 %in% c("SP4210.all", "SP4211.all") |
                                     ind2 %in% c("SP4210.all", "SP4211.all"))
          
          # Extract a subset of 25,000 samples of non-Guam samples
          other_comparisions <- pairwiseDiff %>%
                                filter(!(ind1 %in% c("SP4210.all", "SP4211.all") |
                                         ind2 %in% c("SP4210.all", "SP4211.all"))) %>%
                                filter(noSites > 20000) %>%
                                sample_n(25000)
          
          pairwiseDist_plt <- list(guam_comparisons %>%
                                   mutate(type = ifelse(ind1 == "SP4210.all" & ind2 == "SP4211.all", "intra-Guam", "inter-Guam")),
                                   other_comparisions %>%
                                   mutate(type = "other")) %>%
                              bind_rows() %>%
                              mutate(type = factor(type, levels = c("intra-Guam", "inter-Guam", "other"))) %>%
                              mutate(dummyX = 1) %>%
                              ggplot(aes(x = dummyX, y = fracDiff, colour = type, alpha = type)) +
                              geom_jitter(size = 2.5) +
                              labs(y = "proportion of pairwise differences") +
                              scale_y_continuous(limits = c(0, 0.5),
                                                 labels = scales::percent_format(accuracy = 1)) +
                              scale_alpha_manual(values = c(1, 0.7, 0.25)) +
                              theme_classic(base_size = 11) +
                              theme(legend.position = "top",
                                    axis.title.x = element_blank(),
                                    axis.text.x = element_blank(),
                                    axis.ticks.x = element_blank())
          
          ggsave("{output}", plot = pairwiseDist_plt,
                 width = 160, height = 120, units = "mm", useDingbats = F)
        """)
