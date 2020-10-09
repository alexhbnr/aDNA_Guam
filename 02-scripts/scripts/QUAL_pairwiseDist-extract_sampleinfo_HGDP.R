library(data.table)
library(tidyverse)

# Extract samples present in Reich data set 
samplelist <- fread(snakemake@params[["samplelist"]]) %>%
              pull("Sample name")
reich_dataset <- fread(snakemake@params[["reich_ind"]],
					   header = F, col.names = c("sample", "sex", "population")) %>%
				 rownames_to_column(var = "sampleID") %>%
				 filter(sample %in% samplelist)

# Load individual information
sample_info <- fread(snakemake@params[["pedigree"]]) %>%
               mutate(ind1 = str_c("HGDP", str_pad(ind1, 5, pad = "0")),
                      ind2 = str_c("HGDP", str_pad(ind2, 5, pad = "0")))

# Load pairwise distances
pairwise_diff <- fread(snakemake@input[["diff"]]) %>%
                 head(799) %>%
				 mutate(sample2 = names(.)) %>%
				 pivot_longer(-sample2, names_to = "sample1", values_to = "PWdiff") %>%
				 filter(PWdiff > 0)

# Identify families
pedigree_map <- map_dfr(1:nrow(sample_info), function(i) {
					r <- sample_info[i,]
					if (r$ind1 %in% reich_dataset$sample & r$ind2 %in% reich_dataset$sample) {
					  if (match(r$ind1, reich_dataset$sample) < match(r$ind2, reich_dataset$sample)) {
					    tibble(sample2 = r$ind1,
					  		 sample1 = r$ind2,
					  		 relationship = r$relationship)
					  } else {
					    tibble(sample2 = r$ind2,
					  		 sample1 = r$ind1,
					  		 relationship = r$relationship)
					  }
					} else {
					    tibble(sample2 = r$ind1,
					  		   sample1 = r$ind2,
					  		   relationship = "absent")
					}
				  })

# Add information to pairwise
pairwise_diff_detail <- pairwise_diff %>%
                        left_join(reich_dataset %>%
                        		  select(sample1 = sample, pop1 = population), by = "sample1") %>%
                        left_join(reich_dataset %>%
                        		  select(sample2 = sample, pop2 = population), by = "sample2") %>%
                        mutate(popdiff = ifelse(pop1 == pop2, "intra", "inter")) %>%
                        left_join(pedigree_map, by = c("sample1", "sample2")) %>%
                        mutate(relationship = ifelse(is.na(relationship), "unrelated", relationship))

save(pairwise_diff_detail, reich_dataset, file = snakemake@output[[1]])
