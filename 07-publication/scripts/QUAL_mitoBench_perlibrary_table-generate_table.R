library(data.table)
library(tidyverse)

summary_table <- fread(snakemake@params[["summary"]])

damageprofiles <- fread(snakemake@params[["damage"]]) %>%
                  rename(`id` = `sample`) %>%
                  mutate(`sample` = str_match(id, "(SP421[01])-([A-Z][0-9]+)_([a-z]+)")[,2],
                         `library` = str_match(id, "(SP421[01])-([A-Z][0-9]+)_([a-z]+)")[,3],
                         `data type` = str_match(id, "(SP421[01])-([A-Z][0-9]+)_([a-z]+)")[,4],
                         `data type` = recode(`data type`, `all` = "all reads",
                                                           `deam` = "deaminated only reads"),
                         End = str_replace(End, "p", "'"),
                         Pos = recode(Pos, `1` = "1st", `2` = "2nd", `3` = "3rd"),
                         collabel = str_c(End, " ", Pos, " site [%]"),
                         freq = round(freq * 100, 1)) %>%
                  select(sample, library, `data type`, collabel, freq) %>%
                  pivot_wider(names_from = collabel, values_from = freq)

supptable <- summary_table %>%
             select(sample, library, `data type` = readType, `number of unique reads` = `number of reads`,
                    `mean coverage [fold]` = `mean coverage`, `sites >= 5-fold coverage`,
                    `haplogroup`, `haplogroup quality`, `proportion of contamination DNA (contamMix)`) %>%
             left_join(damageprofiles, by = c("sample", "library", "data type")) %>%
             mutate(`included in per-sample analysis` = as.numeric(str_extract(`proportion of contamination DNA (contamMix)`, "[0-9]+\\.*[0-9]*")) < 25,
                    `included in per-sample analysis` = ifelse(`included in per-sample analysis`, "yes", "no")) %>%
             select(sample:`sites >= 5-fold coverage`, starts_with("5'"), starts_with("3'"),
                    haplogroup:`proportion of contamination DNA (contamMix)`, `included in per-sample analysis`) %>%
             arrange(sample, library, `data type`)

fwrite(supptable, sep = "\t", file = snakemake@output[[1]])
