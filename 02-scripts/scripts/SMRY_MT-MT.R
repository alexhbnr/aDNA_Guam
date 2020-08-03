library(ape)
library(data.table)
library(tidyverse)

fasta_fns <- list.files(snakemake@params[["fastadir"]], full.names = T)
snpad_ns <- map_dfr(fasta_fns, function(fn) tibble(sample = str_replace(basename(fn), ".fa", ""),
                                                  `Ns in snpAD consensus sequence` = sum(read.dna(fn, "fasta", as.character = T) == "n")))

hsd <- map_dfr(snpad_ns$sample, function(s) fread(str_c(snakemake@params[["sampledir"]], "/", s, "/", s, ".hsd"),
                                                  select=c("SampleID", "Found_Polys", "Remaining_Polys"))) %>%
       mutate(Remaining_Polys = str_replace_all(Remaining_Polys, "([0-9]+N |[0-9]+N$)", ""),
              Found_Polys = str_replace_all(Found_Polys, " ", ";"),
              Remaining_Polys = str_replace_all(Remaining_Polys, " ", ";"),
              Remaining_Polys = str_replace(Remaining_Polys, ";$", ""))
       

mt <- fread(snakemake@params[["summary"]]) %>%
      mutate(libraryID = str_match(sample, "SP[0-9]+-([A-Z][0-9]+)_[a-z]+")[, 2],
             `data type` = recode(str_match(sample, "SP[0-9]+-[A-Z][0-9]+_([a-z]+)")[, 2],
                                  `all`  = "all reads", `deam` = "deaminated only reads"),
             `sampleID` = recode(str_match(sample, "(SP[0-9]+)-[A-Z][0-9]+_[a-z]+")[, 2],
                                 `SP4210` = "RBC1", `SP4211` = "RBC2")) %>%
      left_join(snpad_ns, by = "sample") %>%
      left_join(hsd, by = c("sample" = "SampleID")) %>%
      select(`sampleID`, `libraryID`, `data type`, `number of reads`, `mean coverage`,
             `Ns in snpAD consensus sequence`,
             haplogroup, `haplogroup quality`,
             `proportion of contamination DNA (contamMix)`,
             `haplogroup contributions (mixEMT)`,
             `haplogroup-defining polymorphisms` = Found_Polys,
             `additional polymorphisms` = Remaining_Polys) %>%
      mutate(`haplogroup contributions (mixEMT)` = ifelse(str_detect(`haplogroup contributions (mixEMT)`, "NA"),
                                                          "did not converge", `haplogroup contributions (mixEMT)`),
             )

fwrite(mt, sep = "\t", file=snakemake@output[[1]])
