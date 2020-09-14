################################################################################
# Project: Ancient human remains from Guam
# Part: Preparing material for publication
# Step: Supplementary table with mitoBench results per library
#
# Alex Huebner, 14/09/20
################################################################################

workdir: "../snakemake_tmp"

rule all:
    input: 
        "../07-publication/supp_tables/SUM_MT_perlibrary.csv"

rule generate_table:
    output:
        "../07-publication/supp_tables/SUM_MT_perlibrary.csv"
    message: "Generate summary table for the mitoBench summary table per library"
    params: 
        summary = "../05-results/QUAL_mitoBench_perlibrary.csv",
        damage = "../05-results/QUAL_mitoBench_perlibrary_CTfreqs.csv"
    script:
        "scripts/QUAL_mitoBench_perlibrary_table-generate_table.R"
