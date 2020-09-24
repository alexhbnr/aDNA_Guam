################################################################################
# Project: Ancient human remains from Guam
# Part: Quality
# Step: Pairwise distance between Guam samples and other samples
#
# Alex Huebner, 24/03/20
################################################################################

from glob import glob
from itertools import combinations
import os
import re

import numpy as np
import pandas as pd

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/"

rule all:
    input:
        "analysis/qual/pairwiseDist/pairwiseDist_1240K.txt.gz"

#### Based on data set used in publication #####################################

rule calculate_pairwisedist_1240K:
    output:
        "analysis/qual/pairwiseDist/pairwiseDist_1240K.txt.gz"
    message: "Calculate the pairwise distance of Guam samples with respect to other populations in the geographical region"
    params:
        prefix_1240K = "/home/irina_pugach/aDNA_Indonesia_Guam/New_Guam/Data/Merge/Clean.HO-MOD.ancGuamALL.HO-ANC.YRI.FRE"
    run:
        # Read individual file
        inds = pd.read_csv(params.prefix_1240K + ".ind", sep="\t",
                           header=None, names=['individual', 'sex', 'population'])

        # Read Genotype file
        geno = np.genfromtxt(params.prefix_1240K + ".geno", delimiter=1)

        pairwise_diff_results = []
        pairwise_differences = combinations(range(inds['individual'].shape[0]), 2)
        last_pw = inds['individual'][0]
        print(last_pw)
        for pw in pairwise_differences:
            if last_pw != inds['individual'][pw[0]]:
                print(inds['individual'][pw[0]])
                last_pw = inds['individual'][pw[0]]
            nonmissing = geno[:, [pw[0], pw[1]]][np.bitwise_and(geno[:, pw[0]] != 9, geno[:, pw[1]] != 9), :].astype(int)
            for i in range(1):
                nonmissing[nonmissing[:, i] == 1, i] = np.random.choice([0, 2], np.sum(nonmissing[nonmissing[:, i] == 1, i]))
            pairwise_diff_results.append((inds['individual'][pw[0]],
                                          inds['individual'][pw[1]],
                                          np.sum(nonmissing[:, 0] != nonmissing[:, 1]) / nonmissing.shape[0],
                                          nonmissing.shape[0]))

        pd.DataFrame(pairwise_diff_results, columns=['ind1', 'ind2', 'fracDiff', 'noSites']) \
            .to_csv(output[0], sep="\t", index=False,
                    float_format="%.4f", compression="gzip")

#### Based on 1000Genomes ######################################################

rule tgenomes:
    input:
        "results/1000Genomes_pairwiseDiff.RData"

rule download_samplelist:
    output:
        "documentation/1000Genomes_samplelist.xlsx"
    message: "Download sample list of 1000Genomes samples"
    params:
        url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx"
    shell:
        "wget -O {output} {params.url}"

rule calculate_pairwise_differences_tgenomes:
    input:
        "documentation/1000Genomes_samplelist.xlsx",
    output:
        "analysis/qual/pairwise_differences/HGDP_pairwisediff.csv"
    message: "Calculate pairwise differences between all HGDP samples present in Reich dataset"
    shell:
        "julia --threads 16 scripts/scripts/QUAL_pairwiseDist-pairwiseDist.jl"

rule extract_sampleinfo_tgenomes:
    input:
        samplelist = "documentation/1000Genomes_samplelist.xlsx",
        diff = "analysis/qual/pairwise_differences/HGDP_pairwisediff.csv"
    output:
        "results/1000Genomes_pairwiseDiff.RData"
    message: "Extract information about relationship and add to pairwise differences"
    params:
        reich_ind = "/mnt/ancient/ModernHuman/ReichLab/reich_public_geno_v42.4/v42.4.1240K.ind",
        pedigree = "documentation/20130606_g1k.ped"
    script:
        "scripts/QUAL_pairwiseDist-extract_sampleinfo.R"

rule calculate_pairwisedist_ancient:
    output:
        pwdiff = "analysis/qual/pairwise_differences/ancientSamples_pairwiseDiff.csv",
        numsites = "analysis/qual/pairwise_differences/ancientSamples_numsites.csv",
        window_stats = "results/QUAL_pairwiseDist_perwindow.csv"
    message: "Calculate pairwise differences between all ancient samples used in the study"
    params:
        prefix = "/home/irina_pugach/aDNA_Indonesia_Guam/New_Guam/Data/for_Alex/Clean.HO-MOD.ancGuamALL.HO-ANC.YRI.FRE.Liangdao.MEGA",
        popoverview = "documentation/population_overview_SNPdata.csv"
    shell:
        "julia --threads 16 scripts/scripts/QUAL_pairwiseDist-pairwiseDist_ancient.jl"

rule extract_sampleinfo_ancient:
    input:
        pwdiff = "analysis/qual/pairwise_differences/ancientSamples_pairwiseDiff.csv",
        numsites = "analysis/qual/pairwise_differences/ancientSamples_numsites.csv",
        window_stats = "results/QUAL_pairwiseDist_perwindow.csv"
    output:
        "results/ancient_pairwiseDiff.RData"
    message: "Extract information about relationship and add to pairwise differences"
    params:
        prefix = "/home/irina_pugach/aDNA_Indonesia_Guam/New_Guam/Data/for_Alex/Clean.HO-MOD.ancGuamALL.HO-ANC.YRI.FRE.Liangdao.MEGA",
        popoverview = "documentation/population_overview_SNPdata.csv"
    script:
        "scripts/QUAL_pairwiseDist-extract_sampleinfo_ancient.R"
