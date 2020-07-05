library(data.table)
library(tidyverse)
library(patchwork)

# Load shotgun runs
samples <- c("SP4210", "SP4211")
libdirs <- map(samples, function(s) list.files(path = paste0(snakemake@params[["dir"]], s),
                                               pattern = ".+-160818_SN7001204_0542_lane2$",
                                               full.names = T)) %>%
           unlist()

misincorporation <- map_dfr(libdirs, function(fn) {
                      fread(paste(fn, "misincorporation.txt", sep = "/"),
                            select = c(1:8, 10:21)) %>%
                      pivot_longer(`G>A`:`G>T`, names_to = "subst", values_to = "count") %>%
                      pivot_longer(A:T, names_to = "refbase", values_to = "refcount") %>%
                      filter(str_sub(subst, 1, 1) == refbase) %>%
                      group_by(End, Pos, subst) %>%
                      summarise(frac = sum(count) / sum(refcount)) %>%
                      filter(Pos <= 15) %>%
                      mutate(sample = str_match(basename(fn), "(SP[0-9]+)-[A-Z][0-9]+-[0-9]+_.+")[,2])
                    }) %>%
                    select(sample, end = End, pos = Pos, subst, frac) %>%
                    mutate(plot_pos = ifelse(end == "3p", pos * (-1), pos),
                           end = factor(ifelse(end == "3p", "3'", "5'"), levels = c("5'", "3'"))) %>%
                    left_join(tibble(subst = c("C>T", "G>A", "A>C", "A>G", "A>T", "C>A", "C>G", "G>C", "G>T", "T>A", "T>C", "T>G"),
                                     col_cat = c("C>T", "G>A", rep("other", 10))), by = "subst")

plot_damage <- function(misinc, smp, xaxis_title = T) {
  p <- filter(misinc, sample == smp) %>%
       ggplot(aes(x = plot_pos, group = subst)) +
       geom_line(aes(colour = col_cat, y = frac)) +
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
       theme_classic(base_size = 10) +
       theme(legend.position = "top",
             plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
             axis.text.x = element_text(size = 6),
             strip.background.x = element_blank())
  if (!xaxis_title) p + theme(axis.title.x = element_blank())
  else p
}

smiley_plt <- plot_damage(misincorporation, "SP4210", xaxis_title = F) +
              plot_damage(misincorporation, "SP4211", xaxis_title = T) +
              guide_area() +
              plot_layout(ncol = 1, heights = c(1, 1, 0.2), guides = "collect")

ggsave(snakemake@output[[1]], plot = smiley_plt,
       width = 160, height = 100, units = "mm", useDingbats = F)
