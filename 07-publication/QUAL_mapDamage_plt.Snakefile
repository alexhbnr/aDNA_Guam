################################################################################
# Project: Ancient human remains from Guam
# Part: Preparing material for publication
# Step: Plot for aDNA damage patterns
#
# Alex Huebner, 28/02/20
################################################################################

workdir: "../snakemake_tmp"

rule all:
    input: 
        "../07-publication/supp_figures/QUAL_summary_damage.pdf",
        "../07-publication/supp_figures/QUAL_shotgun_damage.pdf"

rule summarise_damageprofiler:
    output:
        "../05-results/QUAL_mapDamage.RData"
    message: "Summarise the damage profile of all libraries for the last 15 position on each end of the read"
    params: 
        dir = "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/analysis/logs/damageprofiler/"
    script:
        "scripts/QUAL_mapDamage_plt-summarise_damageprofiler.R"

rule plot_damage_summary:
    input:
        "../05-results/QUAL_mapDamage.RData"
    output:
        "../07-publication/supp_figures/QUAL_summary_damage.pdf"
    message: "Plot summary of damage plot"
    script:
        "scripts/QUAL_mapDamage_plt-plot_damage_summary.R"

rule plot_damage_shotgun:
    input:
        "../05-results/QUAL_mapDamage.RData"
    output:
        "../07-publication/supp_figures/QUAL_shotgun_damage.pdf"
    message: "Plot summary of damage plot"
    params:
        array = "../01-documentation/capture_arrays.tsv"
    script:
        "scripts/QUAL_mapDamage_plt-plot_damage_shotgun.R"
