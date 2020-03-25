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

rule calculate_pairwisedist_1240K:
    output:
        "analysis/qual/pairwiseDist/pairwiseDist_1240K.txt.gz"
    message: "Calculate the pairwise distance of Guam samples with respect to other populations in the geographical region"
    params:
        prefix_1240K = "/home/irina_pugach/aDNA_Indonesia_Guam/New_Guam/Data/Merge/Clean.HO-MOD.ancGuamALL.HO-ANC.YRI.FRE"
    run:
        # Read individual file
        inds = pd.read_csv(prefix_1240K + ".ind", sep="\t",
                           header=None, names=['individual', 'sex', 'population'])

        # Read Genotype file
        geno = np.genfromtxt(prefix_1240K + ".geno", delimiter=1)

        pairwise_diff_results = []
        pairwise_differences = combinations(range(inds['individual'].shape[0]), 2)
        last_pw = inds['individual'][pw[0]]
        for pw in pairwise_differences:
            if last_pw != inds['individual'][pw[0]]:
                print(inds['individual'][pw[0]])
                last_pw = inds['individual'][pw[0]]
            nonmissing = geno[:,[pw[0], pw[1]]][np.bitwise_and(geno[:, pw[0]] != 9, geno[:, pw[1]] != 9),:]
            pairwise_diff_results.append((inds['individual'][pw[0]],
                                          inds['individual'][pw[1]],
                                          np.sum(nonmissing[:,0] != nonmissing[:,1]) / nonmissing.shape[0],
                                          nonmissing.shape[0]))

        pd.DataFrame(pairwise_diff_results, columns=['ind1', 'ind2', 'fracDiff', 'noSites']) \
            .to_csv(output[0], sep="\t", index=False,
                    float_format="%.4f", compression="gzip")
