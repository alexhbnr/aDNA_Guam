################################################################################
# Project: Ancient human remains from Guam
# Part: Summary
# Step: Summary of the MT results per library
#
# Alex Huebner, 18/11/19
################################################################################

from snakemake.utils import R

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/"

rule MT:
    output:
        "results/SUM_MT.csv"
    message: "Summarise the results from the mitoBench ancient mtDNA pipeline"
    params: 
        summary = "analysis/qual/mtDNA_contamination/summary_table.csv"
    run:
        R("""
        library(data.table)
        library(tidyverse)

        mt <- fread("{params.summary}") %>%
              mutate(`library` = str_match(sample, "SP[0-9]+-([A-Z][0-9]+)_[a-z]+")[, 2],
                     `filter` = str_match(sample, "SP[0-9]+-[A-Z][0-9]+_([a-z]+)")[, 2],
                     `sample` = str_match(sample, "(SP[0-9]+)-[A-Z][0-9]+_[a-z]+")[, 2]) %>%
              select(`sample`, `library`, `filter`, `number of reads`, `mean coverage`,
                     `sites >= 5-fold coverage`, `mode of read length`,
                     `% of Ns in snpAD consensus sequence`,
                     haplogroup, `haplogroup quality`,
                     `proportion of contamination DNA (contamMix)`,
                     `haplogroup contributions (mixEMT)`) %>%
              mutate(`haplogroup contributions (mixEMT)` = ifelse(str_detect(`haplogroup contributions (mixEMT)`, "NA"),
                                                                  "did not converge", `haplogroup contributions (mixEMT)`))
        
        fwrite(mt, sep = "\t", file="{output}")
        """)
