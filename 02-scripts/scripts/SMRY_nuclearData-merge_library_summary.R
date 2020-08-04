library(data.table)
library(tidyverse)

# Information on number of reads
nuclear_capture <- fread(snakemake@params[['libdata']]) %>%
                   filter(array %in% c("390K", "390Ksupp")) %>%
                   group_by(sampleID, libraryID) %>%
                   summarise(across(c(`raw reads`:`unique reads`, `deaminated reads`), sum)) %>%
                   mutate(`% mapped L35MQ25` = round(`reads w/ MQ25` * 100 / `reads w/ L35`, 2),
                          `average_duplications` = round(`reads w/ MQ25` / `unique reads`, 2),
                          sampleID = recode(sampleID, SP4210 = "RBC1", SP4211 = "RBC2"))

# Conditional substitution
load(snakemake@input[["condsubst"]])

# Contamination estimate based on X heterogeneity in males
x_het <- fread(snakemake@params[["xhet"]]) %>%
         filter(filter == "all") %>%
         mutate(cont = str_c(round(`point estimate` * 100, 2), "% (SE ", format(`standard error`, digits=2),")"),
                cont = if_else(is.na(cont), "too little data", cont)) %>%
         select(libraryID = library, `cont. estimate based on chrX heterogeneity in males` = cont)

nuclear_capture %>%
left_join(condsubst_cont %>%
          ungroup() %>%
          select(-sampleID), by = "libraryID") %>%
left_join(x_het, by = "libraryID") %>%
fwrite(., sep = "\t", file = snakemake@output[[1]])
