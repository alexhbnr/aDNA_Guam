library(data.table)
library(tidyverse)
library(readxl)

# Read XLSX data
lab_values <- read_xlsx(snakemake@params[["qpcr"]]) %>%
              select(1, 8:11)
names(lab_values) <- as.character(lab_values[1,])
lab_values <- lab_values[2:nrow(lab_values),]


# Identify Guam samples
samples <- list.files("analysis", full.names = T)
samples <- samples[str_detect(basename(samples), "^SP[0-9]+")]

# Identify list of shotgun sequencing runs
shotgun_runs <- fread(snakemake@params[["capture_array"]]) %>%
                filter(array == "shotgun") %>%
                pull(seqrun_id)

# Identify sequencing files for shotgun runs
seqruns <- sapply(samples, function(s) str_replace_all(list.files(s, pattern = "\\.uniq\\.L35MQ25\\.deam.bam$"),
                                                       "\\.uniq\\.L35MQ25\\.deam.bam$", ""), USE.NAMES = F) %>%
           unlist() %>%
           .[str_split_fixed(., "-", n=3)[,3] %in% shotgun_runs]
          
# Generate sample scaffold
overview_libraries <- tibble(id = seqruns,
                            sampleID = str_match(seqruns, "(SP[0-9]+)-[A-Z][0-9]+-[0-9]+_[A-Z]+[0-9]+_[0-9]+_lane[0-9]")[,2],
                            libraryID = str_match(seqruns, "SP[0-9]+-([A-Z][0-9]+)-[0-9]+_[A-Z]+[0-9]+_[0-9]+_lane[0-9]")[,2],
                            seqrunID = str_match(seqruns, "SP[0-9]+-[A-Z][0-9]+-([0-9]+_[A-Z]+[0-9]+_[0-9]+_lane[0-9])")[,2])

# Summarise analyzeBAM.py
analyzeBAM <- map_df(seqruns, function(s) {
                fread(paste0(snakemake@params[["analyzeBAM_dir"]], s, "-L35MQ25.txt"), skip = 7) %>%
                select(id = sample,
                        `raw reads` = total_reads,
                        `reads w/ L35` = ML_pass,
                        `reads w/ MQ25` = MQ_pass)
            }) %>%
            mutate(`% mapped L35MQ25` = round(`reads w/ MQ25` * 100 / `reads w/ L35`, 3))

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

# Summarise number of reads overlapping with the capture array
## Summarise DamageProfiler
subst_freq <- fread(snakemake@params[["damageprofiler"]],
                    header = T, select = c(2, 3, 22)) %>%
              filter(library %in% seqruns) %>%
              rename(id = `library`,
                     `5'CT` = `1`,
                     `3'CT` = `-1`) %>%
              mutate(`5'CT` = round(`5'CT` * 100, 2),
                     `3'CT` = round(`3'CT` * 100, 2))

## Summarise filterBAM.py
filtered_reads <- fread(snakemake@params[["noreads_deam"]]) %>%
                  filter(id %in% seqruns) %>%
                  rename(`deaminated L35MQ25 reads` = nSamples)

overview_libraries %>%
select(-seqrunID) %>%
left_join(analyzeBAM, by = "id") %>%
left_join(dedup, by = "id") %>%
mutate(`average_duplications` = round(`reads w/ MQ25` / `unique reads`, 2)) %>%
left_join(filtered_reads, by = "id") %>%
left_join(subst_freq, by = "id") %>%
mutate(sampleID = recode(sampleID, SP4210 = "RBC1", SP4211 = "RBC2")) %>%
left_join(lab_values, by = c("libraryID" = "#Index library ID")) %>%
select(sampleID:libraryID,
       `Amount powder (mg)`:`# qPCR (spike-in)`,
       `raw reads`:`deaminated L35MQ25 reads`,
        `5' C>T frequency` = `5'CT`,
        `3' C>T frequency` = `3'CT`) %>%
fwrite(., sep = "\t", file = snakemake@output[[1]])
