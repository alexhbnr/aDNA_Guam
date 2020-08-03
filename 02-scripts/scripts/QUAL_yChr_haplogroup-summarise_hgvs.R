library(data.table)
library(tidyverse)

# Extract information of ISOGG SNPs used by yHaplo
isogg_snps <- fread(str_c(snakemake@params[["dir"]], "/isogg.snps.unique.2016.01.04.txt"),
                    col.names = c("SNP", "haplogroup", "POS", "substitution"), select = c(1:4)) %>%
              as_tibble() %>%
              separate(substitution, c("ancA", "derA"), sep = "->")

# Extract SNPs contributing to path used for haplogroup inference
path <- readLines(str_c(snakemake@params[["dir"]], "/paths.ancGuam.Y.txt")) %>%
        str_split_fixed(., " \\| ", n = 2) %>%
        .[,2] %>%
        str_split(., " ")
path_snps <- map(path[[1]], function(hg) {
               unlist(str_split(str_replace(hg, "[A-Za-z0-9\\-]+:", ""), ","))
             }) %>%
             unlist(.)

yHaplo_res <- fread(str_c(snakemake@params[["dir"]], "/../ancGuam.Y.vcf.gz")) %>%
              left_join(isogg_snps, by = "POS") %>%
              filter(ALT != "." | !is.na(SNP)) %>%
              mutate(obsA = if_else(str_sub(`SP4210.all`, 1, 1) == "0", REF, ALT),
                     type = case_when(
                              is.na(ancA) ~ "not defining",
                              obsA == ancA ~ "ancestral",
                              obsA == derA ~ "derived"
                     ),
                     coverage = str_split_fixed(`SP4210.all`, ":", n=2)[,2],
                     ALT = case_when(
                       ALT == "." & REF == ancA ~ derA,
                       ALT == "." & REF == derA ~ ancA,
                       ALT != "." ~ ALT,
                     ),
                     hgvs = str_c("hg19 chrY.g", POS, REF, ">", ALT),
                     `defining yHaplo path` = if_else(SNP %in% path_snps, "yes", "no")) %>%
              select(hgvs, SNP, haplogroup, `observed allele` = obsA, 
                     `ISOGG allele type` = type, `defining yHaplo path`, coverage)

fwrite(yHaplo_res, sep = "\t", file = snakemake@output[[1]], na = "-")
