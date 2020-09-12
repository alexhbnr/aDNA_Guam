################################################################################
# Project: Ancient human remains from Guam
# Part: Preparation
# Step: Table of capture arrays per sequencing run
#
# Alex Huebner, 17/11/19
################################################################################

import pandas as pd

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam"

rule capture_arrays:
    output:
        "documentation/capture_arrays.tsv"
    message: "Generate list of capture arrays used for each sequencing run"
    run:
        capture_arrays = {'160818_SN7001204_0542_lane2': 'shotgun',
                          '160906_M02279_0022_lane1': 'MT',
                          '161208_SN7001204_0563_lane1': '390Ksupp',
                          '161208_SN7001204_0564_lane1': '390K',
                          '170630_D00829_0053_lane1': 'shotgun',
                          '170822_M02279_0141_lane1': 'MT',
                          '170906_D00829_0070_lane1': 'shotgun',
                          '171020_D00829_0083_lane1': '390K',
                          '171020_D00829_0084_lane1': '390Ksupp',
                          '171110_D00829_0090_lane1': '390K',
                          '171110_D00829_0090_lane2': '390Ksupp',
                          '171110_D00829_0091_lane1': '390Ksupp',
                          '171110_D00829_0091_lane2': '390Ksupp',
                         }
        pd.DataFrame.from_dict(capture_arrays, columns=['array'], orient="index"). \
                to_csv(output[0], sep="\t", index=True, index_label="seqrun_id")
