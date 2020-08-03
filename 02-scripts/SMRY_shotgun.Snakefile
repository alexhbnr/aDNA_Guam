################################################################################
# Project: Ancient human remains from Guam
# Part: Summary
# Step: Summary on major information on the shotgun sequencing data
#
# Alex Huebner, 05/07/20
################################################################################

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/"

rule summary_shotgun:
    output:
        "results/SUM_shotgun.csv"
    message: "Summarise the information on the shotgun sequencing data"
    params:
        capture_array = "documentation/capture_arrays.tsv",
        qpcr = "documentation/overview_Guam_sequencingdata_MMeyer.xlsx",
        analyzeBAM_dir = "analysis/logs/sumstats/",
        dedup_dir = "analysis/logs/dedup/",
        damageprofiler = "analysis/logs/5p3p_substitution_summary.txt",
        noreads_deam = "analysis/logs/deam/noReads_deaminatedonly.txt"
    script:
        "scripts/SMRY_shotgun-summary_shotgun.R"
