library(data.table)
library(tidyverse)
library(patchwork)

load(snakemake@input[[1]])
shotgun_runs <- fread(snakemake@params[["array"]]) %>%
                filter(array == "shotgun") %>%
                pull(seqrun_id)

# Subset shotgun runs
shotgun <- filter(misincorporation, seqrun %in% shotgun_runs) %>%
           mutate(sample = recode(sample, SP4210 = "RBC1",
                                          SP4211 = "RBC2"),
                  plot_pos = ifelse(end == "3p", pos * (-1), pos),
                  end = factor(ifelse(end == "3p", "3'", "5'"), levels = c("5'", "3'")),
                  lib = str_c(sample, " - ", lib)) %>%
           left_join(tibble(subst = c("C>T", "G>A", "A>C", "A>G", "A>T", "C>A", "C>G", "G>C", "G>T", "T>A", "T>C", "T>G"),
                            col_cat = c("C>T", "G>A", rep("other", 10))), by = "subst")

plot_damage_summary <- function(misinc, smp, xaxis_title = T) {
  p <- filter(misinc, lib == smp) %>%
       ggplot(aes(x = plot_pos, y = frac, colour = col_cat, group = subst)) +
       geom_line() +
       facet_wrap(~ end, ncol = 2, scales = "free_x") +
       labs(x = "position at read end",
            y = "frequency",
            colour = "substitution type",
            fill = "substitution type",
            title = smp) +
       scale_x_continuous(breaks = seq(-15, 15, 1)) +
       scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                          limits = c(0, 0.5)) +
       scale_color_manual(values = c("red", "blue", "grey50")) +
       theme_classic(base_size = 8) +
       theme(legend.position = "right",
             plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
             axis.text.x = element_text(size = 4),
             axis.text.y = element_text(size = 6),
             strip.background.x = element_blank())
  if (!xaxis_title) p + theme(axis.title.x = element_blank())
  else p
}

damage_plots <- vector(mode = "list", length = length(unique(shotgun$lib)))
for (i in 1:length(unique(shotgun$lib))) {
    damage_plots[[i]] <- plot_damage_summary(shotgun, unique(shotgun$lib)[[i]], xaxis_title = T)
}
plt <- wrap_plots(damage_plots, ncol = 4, guides = "collect")

ggsave(snakemake@output[[1]], plot = plt,
       width = 250, height = 180, units = "mm", useDingbats = F)
