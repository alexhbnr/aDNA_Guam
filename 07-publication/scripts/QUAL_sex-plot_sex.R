library(data.table)
library(tidyverse)

obs <- fread(snakemake@params[["obs"]]) %>%
       mutate(readType = recode(readType, all = "all fragments",
                                          deam = "deaminated fragments only"),
              sample = recode(sample, SP4210 = "RBC1",
                                      SP4211 = "RBC2"),
              ratio_corr = X / (X + autosomes))

sex_plt <- ggplot(obs,
                  aes(x = sample, y = ratio_corr)) +
           geom_point(size=2.5) +
           facet_wrap(~ readType, nrow = 2) +
           # Add expected values
           scale_y_continuous(limits = c(0.25, 0.55)) +
           geom_hline(yintercept = 0.5, colour = "grey50", lty = 3, size = 1.05) +
           geom_hline(yintercept = 0.33, colour = "grey50", lty = 3, size = 1.05) +
           geom_text(data = data.frame(x = c(-Inf, rep(Inf, 3)),
                                       y = rep(c(0.31, 0.52), 2),
                                       label = c("male", rep("", 2), "female")),
                      aes(x = x, y = y, label = label),
                      size = 4, hjust = rep(c(-0.1, rep(1.1, 3)), 2), inherit.aes = F) +
           labs(x = "sample",
                y = "coverage ratio X / (X + autosomes)") +
           theme_classic(base_size = 10) +
           theme(legend.position = "none",
                 axis.text.x = element_text(size = 9),
                 panel.border = element_rect(color = "grey90", fill = NA, size = 0.9),
                 strip.text.x = element_text(size = 9))

ggsave(snakemake@output[[1]],
       plot = sex_plt,
       dpi = 300,
       unit = "mm",
       width = 80,
       height = 120,
       useDingbats = F)
