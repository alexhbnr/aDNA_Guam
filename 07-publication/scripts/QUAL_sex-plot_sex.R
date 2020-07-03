library(data.table)
library(tidyverse)

obs <- fread(snakemake@params[["obs"]]) %>%
       mutate(flt = recode(flt, all = "all fragments",
                                deam = "deaminated fragments only"))

female_value <- unique(fread(snakemake@params[["exp"]])$female)

sex_plt <- ggplot(obs,
                  aes(x = sample, y = female)) +
           geom_point(size=2.5) +
           facet_wrap(~ flt, nrow = 2) +
           # Add expected values
           scale_y_continuous(limits = c(0.015, 0.045)) +
           geom_hline(yintercept = female_value, colour = "grey50", lty = 3, size = 1.05) +
           geom_hline(yintercept = 0.5 * female_value, colour = "grey50", lty = 3, size = 1.05) +
           #geom_text(data = data.frame(x = c(-Inf, rep(Inf, 3)),
                                       #y = rep(c(0.95, 0.45) * female_value, 2),
                                       #label = c("female", rep("", 2), "male")),
                      #aes(x = x, y = y, label = label),
                      #size = 4, hjust = c(-0.1, rep(1.1, 3)), inherit.aes = F) +
           labs(x = "sample",
                y = "coverage ratio X / (X + autosomes)") +
           theme_classic(base_size = 10) +
           theme(legend.position = "none",
                 axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                 panel.border = element_rect(color = "grey90", fill = NA, size = 0.9),
                 strip.text.x = element_text(size = 9))

ggsave(snakemake@output[[1]],
       plot = sex_plt,
       dpi = 300,
       unit = "mm",
       width = 80,
       height = 120,
       useDingbats = F)
