library(data.table)
library(tidyverse)
library(patchwork)

# Load data
load(snakemake@params[["dist_tgenomes"]])
pwdiff_tgenomes <- pairwise_diff_detail
load(snakemake@params[["dist_ancient"]])
pwdiff_ancient <- pairwise_diff_detail
rm(pairwise_diff_detail)
guam_pwdiff_window <- fread(snakemake@params[["dist_window"]])

# A: Pairwise differences between HGDP samples 
## Subset unrelated to 100,000 observations
set.seed(0)
n <- 100000
pwdiff_tgenomes_subset <- pwdiff_tgenomes %>%
                          mutate(relationship = ifelse(relationship == "unrelated", str_c(relationship, " -\n", popdiff, " pop."), relationship)) %>%
                          group_by(relationship) %>%
                          sample_n(if(n() < n) n() else n) %>%
                          left_join(group_by(., relationship)
                                    %>% count()) %>%
                          mutate(relationship = str_replace(relationship, "Order", "degree"),
                                 relationship = factor(str_c(str_to_lower(relationship), "\n(n=", n, ")"),
                                                       levels = c("first degree\n(n=65)", "second degree\n(n=23)",
                                                                  "third degree\n(n=34)", "unrelated -\nintra pop.\n(n=100000)",
                                                                  "unrelated -\ninter pop.\n(n=100000)")))
## Plot
pwdiff_tgenomes_plt <-  ggplot(pwdiff_tgenomes_subset,
                               aes(x = relationship, y = PWdiff, group = relationship)) +
                        geom_boxplot() +
                        labs(y = "pairwise difference",
                             x = "relationship between pairs of HGDP samples") +
                        scale_y_continuous(limits = c(0.1, 0.45), breaks = seq(0.1, 0.45, 0.05)) +
                        theme_classic(base_size = 9)

# B: Pairwise difference between ancient samples
pwdiff_tgenomes_unrel <- filter(pwdiff_tgenomes_subset, relationship %in% c("first degree\n(n=65)",
                                                                            "unrelated -\nintra pop.\n(n=100000)",
                                                                            "unrelated -\ninter pop.\n(n=100000)")) %>%
                         group_by(relationship) %>%
                         summarise(meanPWdiff = mean(PWdiff)) %>%
                         ungroup() %>%
                         mutate(relationship = factor(c("first degree", "unrel. - intra pop.", "unrel. - inter pop."),
                                                      levels = c("first degree", "unrel. - intra pop.", "unrel. - inter pop.")))

pwdiff_ancient_plt <- pwdiff_ancient %>%
                      filter(nSites >= 5000) %>%
                      mutate(comparison = ifelse(pop_diff == "inter", "inter population", pop1)) %>%
                      left_join(group_by(., comparison) %>% count(), by = "comparison") %>%
                      mutate(comparison = str_c(comparison, "\n(n=", n, ")"),
                             comparison = ifelse(comparison == "French_Polynesia_150BP\n(n=1)", "French_Poly-\nnesia_150BP\n(n=1)", comparison),
                             comparison = factor(comparison, levels = c("anc_Guam\n(n=1)", "French_Poly-\nnesia_150BP\n(n=1)",
                                                                        "anc_Vanuatu\n(n=212)", "Liangdao\n(n=1)", "Vietnam_N\n(n=26)",
                                                                        "Vietnam_BA\n(n=11)", "anc_Laos\n(n=3)", "Thailand_IA\n(n=7)",
                                                                        "Malaysia_H\n(n=1)", "inter population\n(n=1702)"))) %>%
                      ggplot(aes(x = comparison, y = PWdiff, group = comparison)) +
                      geom_boxplot() +
                      geom_hline(data = pwdiff_tgenomes_unrel,
                                 aes(yintercept = meanPWdiff, colour = relationship),
                                 size = 0.5, lty = 3) +
                      geom_vline(xintercept = 9.5, size = 1, linetype = 3, colour = "grey50") +
                      annotate("text", x = 4.5, y = 0.32, hjust = 0.5, size = 2.5, label = "intra") +
                      annotate("text", x = 10, y = 0.32, hjust = 0.5, size = 2.5, label = "inter") +
                      scale_y_continuous(limits = c(0.1, 0.35), breaks = seq(0.1, 0.35, 0.05)) +
                      labs(x = "comparison within and between ancient populations",
                           y = "pairwise differences",
                           colour = "mean pairwise difference in HGDP samples") +
                      theme_classic(base_size = 9) +
                      theme(axis.text.x = element_text(size = 6),
                            legend.position = "top",
                            legend.title = element_text(size = 8))

pairwisediff_plt <- pwdiff_tgenomes_plt +
                    pwdiff_ancient_plt +
                    plot_layout(nrow = 2) +
                    plot_annotation(tag_levels = "A")


ggsave(snakemake@output[[1]],
       plot = pairwisediff_plt,
       height = 140, width = 160, dpi=300, units = "mm")

