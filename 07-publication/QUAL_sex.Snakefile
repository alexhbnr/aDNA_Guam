################################################################################
# Project: Ancient human remains from Guam
# Part: Preparing material for publication
# Step: Summary table of endogenous DNA of shotgun data
#
# Alex Huebner, 30/06/20
################################################################################

rule all:
    input: 
        "supp_figures/QUAL_sex.pdf"

rule plot_sex:
    output:
        "supp_figures/QUAL_sex.pdf"
    message: "Plot the inferred ratio of reads aligning to the X chromosome compared to the autosomes"
    params: 
        obs = "../05-results/QUAL_observed_ratio_X_to_autosomes_shotgun.csv",
    script:
        "scripts/QUAL_sex-plot_sex.R"
