################################################################################
# Project: Ancient human remains from Guam
# Part: Preparing material for publication
# Step: Plot for aDNA damage patterns
#
# Alex Huebner, 28/02/20
################################################################################

from snakemake.utils import R

workdir: "../snakemake_tmp"


rule all:
    input: 
        "../07-publication/supp_figures/QUAL_summary_damage.pdf"

rule summarise_damageprofiler:
    output:
        "../05-results/QUAL_mapDamage.RData"
    message: "Summarise the damage profile of all libraries for the last 15 position on each end of the read"
    params: 
        dir = "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/analysis/logs/damageprofiler/"
    run:
        R("""
          library(data.table)
          library(tidyverse)
          library(dtplyr)
          
          # Load data
          samples <- c("SP4210", "SP4211")
          libdirs <- map(samples, function(s) list.files(path = paste0("{params.dir}", s),
                                                         full.names = T)) %>%
                     unlist()
          
          misincoporation <- map_dfr(libdirs, function(fn) {{
                               fread(paste(fn, "misincorporation.txt", sep = "/"),
                                     select = c(1:8, 10:20)) %>%
                               pivot_longer(`G>A`:`G>C`, names_to = "subst", values_to = "count") %>%
                               pivot_longer(A:T, names_to = "refbase", values_to = "refcount") %>%
                               lazy_dt() %>%
                               filter(str_sub(subst, 1, 1) == refbase) %>%
                               group_by(End, Pos, subst) %>%
                               summarise(frac = sum(count) / sum(refcount)) %>%
                               filter(Pos <= 15) %>%
                               mutate(libname = basename(fn)) %>%
                               as_tibble()
                             }}) %>%
                             mutate(sample = str_match(libname, "(SP[0-9]+)-[A-Z][0-9]+-[0-9]+_.+")[,2],
                                    lib = str_match(libname, "SP[0-9]+-([A-Z][0-9]+)-[0-9]+_.+")[,2]) %>%
                             select(sample, lib, end = End, pos = Pos, subst, frac)
          save(misincoporation, file = "{output}")
        """)

rule plot_damage_summary:
    input:
        "../05-results/QUAL_mapDamage.RData"
    output:
        "../07-publication/supp_figures/QUAL_summary_damage.pdf"
    message: "Plot summary of damage plot"
    run:
        R("""
          library(tidyverse)
          library(patchwork)
          
          # Load data
          load("{input}")
          
          plot_damage <- function(misinc, smp, legend = "none") {{
            misinc %>%
            filter(sample == smp) %>%
            # Summarise across libraries and calculate 95% CI
            group_by(sample, end, pos, subst) %>%
            summarise(meanFrac = mean(frac),
                      sdFrac = sd(frac),
                      seFrac = qnorm(.975) * sdFrac / sqrt(length(frac))) %>%
            ungroup() %>%
            # Annotate for plotting 
            left_join(tibble(subst = c("C>T", "G>A", "A>C", "A>G", "A>T", "C>A", "C>G", "G>C", "T>A", "T>C", "T>G"),
                             col_cat = c("C>T", "G>A", rep("other", 9))), by = "subst") %>%
            mutate(plot_pos = ifelse(end == "3p", pos * (-1), pos),
                   end = factor(ifelse(end == "3p", "3'", "5'"), levels = c("5'", "3'")),
                   seFrac = ifelse(col_cat == "other", 0, seFrac)) %>%
            # Plot
            ggplot(aes(x = plot_pos, group = subst)) +
            geom_line(aes(colour = col_cat, y = meanFrac)) +
            geom_ribbon(aes(ymin = meanFrac - seFrac,
                            ymax = meanFrac + seFrac,
                            fill = col_cat),
                        alpha = 0.5) +
            facet_wrap(~ end, ncol = 2, scales = "free_x") +
            labs(x = "position at read end",
                 y = "mean fraction",
                 colour = "substitution type",
                 fill = "substitution type",
                 title = smp) +
            scale_x_continuous(breaks = seq(-15, 15, 1)) +
            scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                               limits = c(0, 0.5)) +
            scale_color_manual(values = c("red", "blue", "grey50")) +
            scale_fill_manual(values = c("red", "blue", "grey50")) +
            theme_classic(base_size = 10) +
            theme(legend.position = "top",
                  plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
                  axis.text.x = element_text(size = 6),
                  strip.background.x = element_blank())
          }}
          
          plt <- plot_damage(misincoporation, "SP4210") +
                 plot_damage(misincoporation, "SP4211", legend = "right") +
                 guide_area() +
                 plot_layout(ncol = 1, heights = c(1, 1, 0.2), guides = "collect")
          
          ggsave("{output}",
                 height = 120, width = 160, units = "mm", useDingbats = F)
        """)
