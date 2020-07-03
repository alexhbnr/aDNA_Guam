library(tidyverse)
library(patchwork)

# Load data
load(snakemake@input[[1]])

plot_damage <- function(misinc, smp, legend = "none") {
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
}

smiley_plt <- plot_damage(misincoporation, "SP4210") +
              plot_damage(misincoporation, "SP4211", legend = "right") +
              guide_area() +
              plot_layout(ncol = 1, heights = c(1, 1, 0.2), guides = "collect")

ggsave(snakemake@output[[1]], plot = smiley_plt,
       width = 160, height = 100, units = "mm", useDingbats = F)

