library(data.table)
library(tidyverse)

condsubst <- fread(snakemake@params[["freq"]],
                   header = T, select = c(1, 2, 3, 6),
                   col.names = c("sample", "seqrun", "cond5p", "cond3p")) %>%
             mutate(lib = str_match(seqrun, "SP421[01]-([A-Z][0-9]+)-[0-9]+_[A-Z]+[0-9]+_[0-9]+_lane[0-9]")[,2]) %>%
             pivot_longer(`cond5p`:`cond3p`, names_to = "end", values_to = "freq") %>%
             group_by(sample, lib, end) %>%
             summarise(meanFreq = mean(freq),
                       stdFreq = sd(freq))

allsubst <- fread(snakemake@params[["subst"]],
                    header = T,
                    select = c(1, 2, 3, 22),
                    col.names = c("sample", "seqrun", "all5p", "all3p")) %>%
            mutate(lib = str_match(seqrun, "SP421[01]-([A-Z][0-9]+)-[0-9]+_[A-Z]+[0-9]+_[0-9]+_lane[0-9]")[,2]) %>%
            pivot_longer(`all5p`:`all3p`, names_to = "end", values_to = "freq") %>%
            group_by(sample, lib, end) %>%
            summarise(meanFreq = mean(freq),
                      stdFreq = sd(freq))

condsubst_cont <- left_join(allsubst %>%
                            select(-stdFreq) %>%
                            pivot_wider(names_from = "end", values_from = "meanFreq"),
                            condsubst %>%
                            select(-stdFreq) %>%
                            pivot_wider(names_from = "end", values_from = "meanFreq"),
                            by = c("sample" = "sample",
                                   "lib" = "lib")) %>%
                  mutate(cont3p = 1 - (all3p / cond3p),
                         cont5p = 1 - (all5p / cond5p)) %>%
                  pivot_longer(cont3p:cont5p, names_to = "readend", values_to = "contamination") %>%
                  mutate(contamination = if_else(contamination < 0, 0, contamination)) %>%
                  group_by(sample, lib, all3p, all5p, cond3p, cond5p) %>%
                  summarise(contamination = mean(contamination)) %>%
                  mutate(sample = recode(sample, SP4210 = "RBC1", SP4211 = "RBC2")) %>%
                  rename(sampleID = sample, libraryID = lib,
                         `C>T frequency 3' end - all reads` = all3p,
                         `C>T frequency 5' end - all reads` = all5p,
                         `C>T frequency 3' end - conditional` = cond3p,
                         `C>T frequency 5' end - conditional` = cond5p,
                         `C>T frequency 3' end - conditional` = cond3p,
                         `C>T frequency 5' end - conditional` = cond5p,
                         `contamination estimate based on deamination` = contamination)
save(condsubst_cont, file = snakemake@output[[1]])

