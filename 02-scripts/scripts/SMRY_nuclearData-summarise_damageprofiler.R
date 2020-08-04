library(data.table)
library(tidyverse)
library(dtplyr)

samples <- list.files(snakemake@params[["dir"]], full.names = T)
samples <- list.files("tmp/damageprofiler", full.names = T)

misincorporation <- map_dfr(samples, function(fn) {
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
                    }) %>%
                    mutate(lib = str_match(libname, "([A-Z][0-9]+)\\.([a-z+])")[,2],
                           type = str_match(libname, "([A-Z][0-9]+)\\.([a-z+])")[,3]) %>%
                    select(lib, type, end = End, pos = Pos, subst, frac)

save(misincorporation, file = snakemake@output[[1]])
