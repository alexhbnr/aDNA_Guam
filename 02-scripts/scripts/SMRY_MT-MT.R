library(ape)
library(data.table)
library(tidyverse)

samples <- list.files(snakemake@params[["sampledir"]], pattern = "^SP")

# Read FastA files
snpad_ns <- map_dfr(samples, function(s) tibble(sample = s,
                                                `Ns in snpAD consensus sequence` = sum(read.dna(str_c(snakemake@params[["sampledir"]], "/", s, "/", s, ".fa"),
                                                                                                "fasta", as.character = T) == "n"))) %>%
            mutate(`data type` = recode(str_match(sample, "(SP421[01])_([a-z]+)")[,3],
                                        `all` = "all reads", `deam` = "deaminated only reads"),
                   `sample` = str_match(sample, "(SP421[01])_([a-z]+)")[,2])
# Read HaploGrep output
hsd <- map_dfr(samples, function(s) fread(str_c(snakemake@params[["sampledir"]], "/", s, "/", s, ".hsd"),
                                          select=c("SampleID", "Found_Polys", "Remaining_Polys"), nrows=1)) %>%
       mutate(Remaining_Polys = str_replace_all(Remaining_Polys, "([0-9]+N |[0-9]+N$)", ""),
              Found_Polys = str_replace_all(Found_Polys, " ", ";"),
              Remaining_Polys = str_replace_all(Remaining_Polys, " ", ";"),
              Remaining_Polys = str_replace(Remaining_Polys, ";$", ""),
              `data type` = recode(str_match(SampleID, "(SP421[01])_([a-z]+)")[,3],
                                   `all` = "all reads", `deam` = "deaminated only reads"),
              `sample` = str_match(SampleID, "(SP421[01])_([a-z]+)")[,2])
       

mt <- fread(snakemake@params[["summary"]]) %>%
      rename(`data type` = readType) %>%
      left_join(snpad_ns, by = c("sample", "data type")) %>%
      left_join(hsd, by = c("sample", "data type")) %>%
      select(`sample`, `data type`,
             `number of unique reads` = `number of reads`,
             `mean coverage [fold]` = `mean coverage`,
             `Ns in snpAD consensus sequence`,
             haplogroup, `haplogroup quality`,
             `proportion of contamination DNA (contamMix)`,
             `haplogroup-defining polymorphisms` = Found_Polys,
             `additional polymorphisms` = Remaining_Polys)

fwrite(mt, sep = "\t", file=snakemake@output[[1]])
