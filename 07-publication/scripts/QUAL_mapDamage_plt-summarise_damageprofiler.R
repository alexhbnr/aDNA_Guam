library(data.table)
library(tidyverse)
library(dtplyr)

# Load data
samples <- c("SP4210", "SP4211")
libdirs <- map(samples, function(s) list.files(path = paste0(snakemake@params[["dir"]], s),
                                               full.names = T)) %>%
           unlist()

misincoporation <- map_dfr(libdirs, function(fn) {
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
                   mutate(sample = str_match(libname, "(SP[0-9]+)-[A-Z][0-9]+-[0-9]+_.+")[,2],
                          lib = str_match(libname, "SP[0-9]+-([A-Z][0-9]+)-[0-9]+_.+")[,2]) %>%
                   select(sample, lib, end = End, pos = Pos, subst, frac)
save(misincoporation, file = snakemake@output[[1]])

