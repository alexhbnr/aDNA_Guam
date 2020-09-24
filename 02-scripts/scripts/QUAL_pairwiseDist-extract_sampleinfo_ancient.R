library(data.table)
library(tidyverse)

# A: Pairwise differences between ancient individuals from same region
## Load information on individuals
inds <- fread(str_c(snakemake@params[["prefix"]], ".ind"),
              col.names = c("ind", "sex", "pop"), header=F)
## Load information on populations
pops <- fread(snakemake@params[["popoverview"]])

## Load pairwise distances
pairwise_diff <- fread(snakemake@input[["pwdiff"]]) %>%
				 mutate(sample2 = colnames(.)) %>%
				 pivot_longer(-sample2, names_to = "sample1", values_to = "PWdiff") %>%
				 filter(PWdiff > 0)

## Load number of sites
numsites <- fread(snakemake@input[["numsites"]]) %>%
				  mutate(sample2 = colnames(.)) %>%
				  pivot_longer(-sample2, names_to = "sample1", values_to = "nSites") %>%
				  filter(nSites > 0)

## Combine data
pairwise_diff_detail <- pairwise_diff %>%
                        left_join(numsites, by = c("sample2", "sample1")) %>%
                        left_join(inds %>%
                                  select(sample2 = ind, pop2 = pop), by = "sample2") %>%
                        left_join(inds %>%
                                  select(sample1 = ind, pop1 = pop), by = "sample1") %>%
                        mutate(pop_diff = if_else(pop2 == pop1, "intra", "inter"))

## Calculate summary statistics on the intra-population PW difference
intra_pwdiff <- filter(pairwise_diff_detail, pop_diff == "intra") %>%
                group_by(pop1) %>%
                summarise(meanPWD = mean(PWdiff),
                          sdPWD = sd(PWdiff),
                          n = length(PWdiff))

## Calculate summary statistics on the inter-population PW difference
inter_pwdiff <- filter(pairwise_diff_detail, pop_diff == "inter") %>%
                summarise(meanPWD = mean(PWdiff),
                          sdPWD = sd(PWdiff),
                          n = length(PWdiff))

save(pairwise_diff_detail, file = snakemake@output[[1]])
