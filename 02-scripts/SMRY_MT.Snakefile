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
        summary = "analysis/qual/mtDNA_contamination/summary_table.csv",
        fastadir = "analysis/qual/mtDNA_contamination/fasta",
        sampledir = "analysis/qual/mtDNA_contamination/sample_stats",
    script:
        "scripts/SMRY_MT-MT.R"
