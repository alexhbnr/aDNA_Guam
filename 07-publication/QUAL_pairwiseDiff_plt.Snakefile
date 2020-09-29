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
        "supp_figures/QUAL_pairwiseDist.pdf"

rule pairwiseDiff_plt:
    output:
        "supp_figures/QUAL_pairwiseDist.pdf"
    message: "Plot the pairwise difference analysis results for the Guam samples"
    params: 
        dist = "../04-analysis/qual/pairwiseDist/pairwiseDist_1240K.txt.gz",
        inds = "../01-documentation/mergedDataset_Clean.HO-MOD.ancGuamALL.HO-ANC.YRI.FRE.ind"
    run:
        R("""
          library(data.table)
          library(tidyverse)
          
          pairwiseDiff <- fread("{params.dist}")
          inds <- fread("{params.inds}",
                        header = F, col.names = c("individual", "sex", "population"))

          pairwiseDiff %<>%
          # Add population info for ind1
          left_join(inds %>%
                    select(ind1 = individual, pop1 = population), by = "ind1") %>%
          # Add population info for ind2
          left_join(inds %>%
                    select(ind2 = individual, pop2 = population), by = "ind2") %>%
          # Infer comparison type
          mutate(type = if_else(pop1 == pop2, "intra", "inter"),
                 type = if_else(pop1 == "anc_Guam" | pop2 == "anc_Guam", paste(type, "Guam", sep="-"),
                                                                         paste(type, "other", sep="-")))

          pairwiseDiff_subset <- pairwiseDiff %>%
                                 group_by(type, guam = str_detect(type, "Guam")) %>%
                                 nest() %>%
                                 mutate(subset = map(data, ~ {if (nrow(.x) > 25000) {
                                                                sample_n(.x, 25000) 
                                                              } else {
                                                                .x
                                                              }})) %>%
                                 select(-data) %>%
                                 unnest(cols = "subset") %>%
                                 ungroup() %>%
                                 select(-guam)

          pairwiseDist_plt <- pairwiseDiff_subset %>%
                              mutate(type = factor(type, levels = c("intra-Guam", "intra-other",
                                                                    "inter-Guam", "inter-other"))) %>%
                              ggplot(aes(x = type, y = fracDiff, colour = type)) +
                              geom_boxplot() +
                              geom_vline(xintercept = 2.5, colour = "black", lty = 3) +
                              annotate("text", x = 1.5, y = 0.5, label = "intra-population") +
                              annotate("text", x = 3.5, y = 0.5, label = "inter-population") +
                              labs(x = "comparison types",
                                   y = "proportion of pairwise different sites") +
                              scale_x_discrete(labels = rep(c("Guam", "other populations"), 2)) +
                              scale_y_continuous(limits = c(0, 0.5),
                                                 labels = scales::percent_format(accuracy = 1)) +
                              scale_colour_manual(values = c("red", rep("black", 3))) +
                              theme_classic() +
                              theme(legend.position = "none")
          
          ggsave("{output}", plot = pairwiseDist_plt,
                 width = 160, height = 120, units = "mm", useDingbats = F)
        """)


rule pairwiseDiff_tgenomes:
    output:
        "../07-publication/supp_figures/QUAL_pairwiseDist_1000Genomes.pdf"
    message: "Plot the pairwise difference analysis results for the Guam samples using the 1000Genomes dataset for comparison"
    params: 
        dist_ancient = "../05-results/ancient_pairwiseDiff.RData",
        dist_tgenomes = "../05-results/1000Genomes_pairwiseDiff.RData",
        dist_window = "../05-results/QUAL_pairwiseDist_perwindow.csv"
    script:
        "scripts/QUAL_pairwiseDist_plt-pairwiseDiff_tgenomes.R"
