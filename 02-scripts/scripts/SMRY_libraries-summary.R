library(data.table)
library(tidyverse)

# Identify Guam samples
samples <- list.files("analysis", full.names = T)
samples <- samples[str_detect(basename(samples), "^SP[0-9]+")]
# Identify sequencing files
seqruns <- sapply(samples, function(s) str_replace_all(list.files(s, pattern = "\\.uniq\\.L35MQ25\\.deam.bam$"),
                                                       "\\.uniq\\.L35MQ25\\.deam.bam$", ""), USE.NAMES = F) %>%
           unlist()

# Generate sample scaffold
overview_libraries <- tibble(id = seqruns,
                             sampleID = str_match(seqruns, "(SP[0-9]+)-[A-Z][0-9]+-[0-9]+_[A-Z]+[0-9]+_[0-9]+_lane[0-9]")[,2],
                             libraryID = str_match(seqruns, "SP[0-9]+-([A-Z][0-9]+)-[0-9]+_[A-Z]+[0-9]+_[0-9]+_lane[0-9]")[,2],
                             seqrunID = str_match(seqruns, "SP[0-9]+-[A-Z][0-9]+-([0-9]+_[A-Z]+[0-9]+_[0-9]+_lane[0-9])")[,2]) %>%
                      # Add array information
                      left_join(fread(snakemake@params[["array"]]),
                                by=c("seqrunID"="seqrun_id"))
# Summarise analyzeBAM.py
analyzeBAM <- map_df(seqruns, function(s) {
                fread(paste0(snakemake@params[["analyzeBAM_dir"]], s, "-L35MQ25.txt"), skip = 7) %>%
                select(id = sample,
                       `raw reads` = total_reads,
                       `reads w/ L35` = ML_pass,
                       `reads w/ MQ25` = MQ_pass)
              })
overview_libraries <- left_join(overview_libraries,
                                analyzeBAM, by = "id")

# Summarise dedup
dedup <- map_df(seqruns, function(s) {
           # Parse dedup logfile
           loginfo <- readLines(paste0(snakemake@params[["dedup_dir"]], s, ".L35MQ25.log"))
           totalreads <- as.numeric(str_match(loginfo[2], "Total reads: ([0-9]+)")[, 2])
           removedreads <- as.numeric(str_match(loginfo[6], "Total removed: ([0-9]+)")[, 2])
           duprate <- as.numeric(str_match(loginfo[7], "Duplication Rate: ([0-9\\.]+)")[, 2])
           tibble(id = s,
                  `unique reads` = totalreads - removedreads)
         })
overview_libraries <- left_join(overview_libraries,
                                dedup, by = "id") %>%
                      mutate(`duplication rate` = round(`reads w/ MQ25` / `unique reads`, 2))
# Summarise number of reads overlapping with the capture array
tf_overlap <- map_df(list.files(snakemake@params[["genotypes_array"]],
                     pattern = "\\.nsnps\\.txt", full.names = T), function(fn) {
                       fread(fn, header = F,
                             col.names = c("id", "array", "nSNPs"))
                    }) %>%
              mutate(experiment = ifelse(str_sub(id, nchar(id) - 3, nchar(id)) == "deam", "deam", "all"),
                     id = str_replace(id, "\\.uniq\\.L35MQ25.*", ""),
                     experiment_array = paste0(array, " SNPs covered by " , experiment)) %>%
              select(-array, -experiment) %>%
              spread(experiment_array, nSNPs, fill = 0)
overview_libraries <- left_join(overview_libraries,
                                tf_overlap %>%
                                select(id, `390K SNPs covered by all`, `390Ksupp SNPs covered by all`, `archaicAdmix SNPs covered by all`),
                                by = "id")

# Summarise DamageProfiler
subst_freq <- fread(snakemake@params[["damageprofiler"]],
                    header = T, select = c(2, 3, 22)) %>%
              rename(id = `library`,
                     `5'CT` = `1`,
                     `3'CT` = `-1`) %>%
              mutate(`5'CT` = round(`5'CT` * 100, 2),
                     `3'CT` = round(`3'CT` * 100, 2))
overview_libraries <- left_join(overview_libraries,
                                subst_freq, by = "id")
average_length <- fread(snakemake@params[["averagelength"]], select = c(2, 5)) %>%
                  rename(id = `library`,
                         `mode read length` = mode)
overview_libraries <- left_join(overview_libraries,
                                average_length, by = "id")

# Summarise filterBAM.py
filtered_reads <- fread(snakemake@params[["noreads_deam"]]) %>%
                  rename(`deaminated reads` = nSamples)
overview_libraries <- left_join(overview_libraries,
                              filtered_reads, by = "id")
overview_libraries <- left_join(overview_libraries,
                                tf_overlap %>%
                                select(id, `390K SNPs covered by deam`, `390Ksupp SNPs covered by deam`, `archaicAdmix SNPs covered by deam`),
                                by = "id")

# Summarise conditional substitution contamination estimate
condsubst_cont <- fread(snakemake@params[["condsubst_cont"]],
                        select = c(2, 3)) %>%
                  rename(id = `library`) %>%
                  mutate(`contamination [%]` = ifelse(contamination > 0, round(contamination * 100, 1), 0)) %>%
                  select(-contamination)
overview_libraries <- left_join(overview_libraries,
                                condsubst_cont, by = "id") %>%
                      select(-id)

# Summarise X heterogenity for males
xheterogeneity <- fread(snakemake@params[["xheterogeneity"]]) %>%
                 mutate(array = "390K",
                        `point estimate` = round(`point estimate` * 100, 2)) %>%
                 filter(filter == "all") %>%
                 select(-`standard error`, -filter) %>%
                 rename(`cont. based on X heterogeneity [%] - all reads` = `point estimate`)
overview_libraries <- left_join(overview_libraries,
                                xheterogeneity,
                                by = c("sampleID" = "sample",
                                       "libraryID" = "library",
                                       "array" = "array")) %>%
                      distinct() %>%
                      arrange(sampleID, libraryID, array, seqrunID)

fwrite(overview_libraries, sep = "\t", file = snakemake@output[[1]])
