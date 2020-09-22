library(data.table)
library(tidyverse)
library(readxl)

# Load 1000Genomes sample info
sample_info <- read_xlsx(snakemake@input[["samplelist"]])

# Extract samples present in Reich data set 
reich_dataset <- fread(snakemake@params[["reich_ind"]],
					   header = F, col.names = c("sample", "sex", "population")) %>%
				 rownames_to_column(var = "sampleID") %>%
				 mutate(hgdpid = str_extract(sample, "[A-Z]+[0-9]+")) %>%
				 filter(hgdpid %in% sample_info$Sample)
hgdp_reich <- sample_info %>%
			  filter(Sample %in% str_extract(reich_dataset$sample, "[A-Z]+[0-9]+"))

# Load pairwise distances
pairwise_diff <- fread(snakemake@input[["diff"]]) %>%
				 mutate(sample2 = colnames(.)) %>%
				 pivot_longer(-sample2, names_to = "sample1", values_to = "PWdiff") %>%
				 filter(PWdiff > 0)

# Identify families
pedigree <- fread(snakemake@params[["pedigree"]]) %>%
			filter(`Individual ID` %in% hgdp_reich$Sample) %>%
			pivot_longer(c(`Paternal ID`:`Maternal ID`, Siblings:`Third Order`), names_to = "rel_type", values_to = "rel_info") %>%
			filter(rel_info != "0")

pedigree_map <- map_dfr(pedigree$`Individual ID`, function(id) {
				  relationships <- filter(pedigree, `Individual ID` == id) %>%
				                   mutate(rel_type = recode(rel_type,
				                 						   `Maternal ID` = "First Order",
				                 						   `Paternal ID` = "First Order",
														   Siblings = "First Order")) %>%
								   separate_rows(rel_info, convert=T)
				  map_dfr(1:nrow(relationships), function(i) {
					r <- relationships[i,]
					if (r$rel_info %in% reich_dataset$hgdpid) {
					  if (match(r$`Individual ID`, reich_dataset$hgdpid) < match(r$rel_info, reich_dataset$hgdpid)) {
					    tibble(sample1 = r$`Individual ID`,
					  		 sample2 = r$rel_info,
					  		 relationship = r$rel_type)
					  } else {
					    tibble(sample1 = r$rel_info,
					  		 sample2 = r$`Individual ID`,
					  		 relationship = r$rel_type)
					  }
					} else {
					    tibble(sample1 = r$`Individual ID`,
					  		 sample2 = r$rel_info,
					  		 relationship = "absent")
					}
				  })
				})

# Add information to pairwise
pairwise_diff_detail <- pairwise_diff %>%
                        left_join(hgdp_reich %>%
                        		  select(sample1 = Sample, pop1 = Population), by = "sample1") %>%
                        left_join(hgdp_reich %>%
                        		  select(sample2 = Sample, pop2 = Population), by = "sample2") %>%
                        mutate(popdiff = ifelse(pop1 == pop2, "intra", "inter")) %>%
                        left_join(pedigree_map, by = c("sample1", "sample2")) %>%
                        mutate(relationship = ifelse(is.na(relationship), "unrelated", relationship))

save(pairwise_diff_detail, hgdp_reich, file = snakemake@output[[1]])
