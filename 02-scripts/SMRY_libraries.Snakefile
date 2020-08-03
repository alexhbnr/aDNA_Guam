################################################################################
# Project: Ancient human remains from Guam
# Part: Summary
# Step: Summary on major information on the libraries
#
# Alex Huebner, 18/11/19
################################################################################

from snakemake.utils import R

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/"

rule summary:
    output:
        "results/SUM_libraries.csv"
    message: "Summarises all information regarding the libraries of the Guam samples"
    params: 
        array = "documentation/capture_arrays.tsv",
        analyzeBAM_dir = "analysis/logs/sumstats/",
        dedup_dir = "analysis/logs/dedup/",
        genotypes_array = "analysis/logs/nSNPs",
        damageprofiler = "analysis/logs/5p3p_substitution_summary.txt",
        averagelength = "analysis/logs/average_length.txt",
        noreads_deam = "analysis/logs/deam/noReads_deaminatedonly.txt",
        condsubst_cont = "analysis/qual/conditionalSubsts_contaminationEstimate.txt",
        xheterogeneity = "analysis/qual/Xheterogeneity/X_contamination.txt"
    script:
        "scripts/SMRY_libraries-summary.R"
