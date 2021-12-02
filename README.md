# Parts of the analyses of "Ancient DNA from Guam and the peopling of the Pacific"

This repository contains the scripts for processing of the sequencing data and some analyses of the
paper

> I Pugach, A Hübner, H Hung, M Meyer, MT Carson, M Stoneking: Ancient DNA from Guam and the peopling of the Pacific. *PNAS* **118** (1) e2022112118 (2021)

DOI: [10.1073/pnas.2022112118](https://doi.org/10.1073/pnas.2022112118)

## Table of Contents

This repository contains all analyses performed by myself, which resulted in either a Supplementary
table or a Supplementary figure.

### Code for Supplementary Tables

#### Table S1: Overview of the shotgun sequencing results

  - [code for generating table](02-scripts/SMRY_shotgun.Snakefile)
  - [table](07-publication/supp_tables/SUM_shotgun.xlsx)

#### Table S2: Overview of the mtDNA-enriched sequencing data

  - [code for generating table](07-publication/QUAL_mitoBench_perlibrary_table.Snakefile)
  - [table](07-publication/supp_tables/SUM_MT_perlibrary.xlsx)

#### Table S3: Overview of the mtDNA analysis

  - [code for generating table](02-scripts/MT_mitoBench.Snakefile)
  - [table](07-publication/supp_tables/SUM_MT.xlsx)

#### Table S4: Overview of the Y-chromosome analysis

  - [code for generating table](02-scripts/QUAL_yChr_haplogroup.Snakefile)
  - [table](07-publication/supp_tables/YCHR_hgvs.xlsx)

#### Table S5: Overview of the 1240K-enriched sequencing data

  - [code for generating table](02-scripts/SMRY_nuclearData.Snakefile)
  - [table](07-publication/supp_tables/SUM_nuclearcapture.xlsx)

#### Table S6: Overview of the number of nuclear SNPs per sample

  - [code for generating table](02-scripts/SMRY_samples.Snakefile)
  - [table](07-publication/supp_tables/SUM_SNPs.xlsx)

### Code for Supplementary Figures

#### Figure S3: Ancient DNA damage in the shotgun sequencing data per library

  - [code for generating data](02-scripts/PREP_filterSamples.Snakefile)
  - [results](05-results/QUAL_mapDamage.RData)
  - [code for plotting](07-publication/QUAL_mapDamage_plt.Snakefile)
  - [figure](07-publication/supp_figures/QUAL_shotgun_damage.pdf)

#### Figure S4: Genetic sex determination from shotgun sequencing data

  - [code for generating data](02-scripts/QUAL_sex_shotgun.Snakefile)
  - [results](05-results/QUAL_observed_ratio_X_to_autosomes_shotgun.csv)
  - [code for plotting](07-publication/QUAL_sex.Snakefile)
  - [figure](07-publication/supp_figures/QUAL_sex.pdf)

#### Figure S5: Comparison of the fraction of pairwise differences among present-day and ancient populations

  - [code for generating data](02-scripts/QUAL_pairwiseDist.Snakefile)
  - [code for plotting](07-publication/QUAL_pairwiseDiff_plt.Snakefile)
  - [figure](07-publication/supp_figures/QUAL_pairwiseDist_1000Genomes.pdf)
