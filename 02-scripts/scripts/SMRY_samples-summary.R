library(data.table)
library(tidyverse)

snp_fns <- list.files(snakemake@params[["dir"]], full.names = T)
print(snp_fns)
snp_summary <- map_df(snp_fns, function(fn) {
                 fread(fn) %>%
                 mutate(array = str_replace(basename(fn), "\\.txt", ""))
               }) %>%
               spread(array, nSNPs) %>%
               mutate(`1240K` = `390K` + `390Ksupplement`,
                      filter = ifelse(filter == "deam", "deaminated only", "all")) %>%
               select(sample, reads = filter,
                      `390K`, `390Ksupplement`, `1240K`, `archaicadmixture`,
                      HumanOrigins, `Affymetrix6.0` = `Affymetrix6`,
                      `Axiom GWAS CEU1` = Axiom_GW_CEU1, MEGA)
fwrite(snp_summary, sep = "\t", file = snakemake@output[[1]])
