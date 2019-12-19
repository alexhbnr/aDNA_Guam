################################################################################
# Project: Ancient human remains from Guam
# Part: Data preparation
# Step: De-multiplex sequencing runs based on perfect-index matching
#
# Alex Huebner, 17/11/19
################################################################################

import os.path

import pandas as pd

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam"

#### Sequencing runs ###########################################################
SEQRUNS = {'160906_M02279_0022_lane1': '/mnt/ngs_data/160906_M02279_0022_000000000-AUCAV_BN_D2869/Bustard/BWA/proc1/s_1_sequence_ancient_human_MT.bam',
           '161208_SN7001204_0563_lane1': '/mnt/ngs_data/161208_SN7001204_0563_AH3WJVBCXY_R_PEdi_D4262_D4264/Bustard/BWA/proc1/s_1_sequence_ancient_hg19_evan.bam',
           '161208_SN7001204_0564_lane1': '/mnt/ngs_data/161208_SN7001204_0564_BH3WKLBCXY_R_PEdi_D4267_D4268/Bustard/BWA/proc1/s_1_sequence_ancient_hg19_evan.bam',
           '170906_D00829_0070_lane1': '/mnt/ngs_data/170906_D00829_0070_AHNVFNBCXY_R_PEdi_F8888_G1774/Bustard/BWA/proc1/s_1_sequence_ancient_hg19_evan.bam',
           '171020_D00829_0083_lane1': '/mnt/ngs_data/171020_D00829_0083_AHWTHYBCXY_R_PEdi_F5903_F5897/Bustard/BWA/proc1/s_1_sequence_ancient_hg19_evan.bam',
           '171020_D00829_0084_lane1': '/mnt/ngs_data/171020_D00829_0084_BHWVGYBCXY_R_PEdi_F5904_G3070/Bustard/BWA/proc1/s_1_sequence_ancient_hg19_evan.bam',
           '171110_D00829_0090_lane1': '/mnt/ngs_data/171110_D00829_0090_AH2GCVBCX2_R_PEdi_F5903_F5904/Bustard/BWA/proc2/s_1_sequence_ancient_hg19_evan.bam',
           '171110_D00829_0090_lane2': '/mnt/ngs_data/171110_D00829_0090_AH2GCVBCX2_R_PEdi_F5903_F5904/Bustard/BWA/proc1/s_2_sequence_ancient_hg19_evan.bam',
           '171110_D00829_0091_lane1': '/mnt/ngs_data/171110_D00829_0091_BH2M2TBCX2_R_PEdi_F5904/Bustard/BWA/proc1/s_1_sequence_ancient_hg19_evan.bam',
           '171110_D00829_0091_lane2': '/mnt/ngs_data/171110_D00829_0091_BH2M2TBCX2_R_PEdi_F5904/Bustard/BWA/proc1/s_2_sequence_ancient_hg19_evan.bam',
           }
################################################################################

#### Libraries #################################################################
LIBRARIES = {'SP4210': ['D5864', 'F7638', 'F7639', 'F7640', 'F7641', 'F8851', 'F8852', 'F8853', 'F8854', 'F8855', 'F8856', 'F8857', 'F8858', 'F5861'],  # SP4210 - RBC1, Ritidian Beach Cave, Guam
             'SP4211': ['D5865', 'F7642', 'F7643', 'F7644', 'F7645', 'F8862', 'F8863', 'F8864', 'F8865', 'F8866', 'F8867'],  # SP4211 - RBC2, Ritidian Beach Cave, Guam
            }
# Generate list of all read groups in BAM files that are to be expected
LIBLIST = [lib for _, v in LIBRARIES.items() for lib in v]
################################################################################

localrules: index_list

rule demultiplexing:
    input: expand("data/splitBAM/{seqrun_id}.done", seqrun_id=SEQRUNS.keys())

rule index_list:
    output:  
        "data/splitBAM/{seqrun_id}.indexlist.txt"
    params: indexlist = lambda wildcards: os.path.dirname(SEQRUNS[wildcards.seqrun_id]) + "/../../../build/lane" + wildcards.seqrun_id[-1] + "/indices.txt"
    run:
        indexlist = pd.read_csv(params.indexlist, sep="\t", index_col=[2])
        # Remove all index combinations that are not from Flores samples
        indexlist = indexlist.loc[LIBLIST].dropna(axis=0, how="all")
        indexlist = indexlist.reset_index()
        indexlist.columns = ['readgroup', 'p7', 'p5']
        indexlist.to_csv(output[0], sep="\t", index=False)

rule sort_BAM_by_name:
    output:
        temp("data/splitBAM/{seqrun_id}.nsorted.bam")
    message: "Sort the BAM file {wildcards.seqrun_id} by name"
    params:
        bam = lambda wildcards: SEQRUNS[wildcards.seqrun_id]
    shell:
        """
        samtools sort \
                -l 0 \
                -n \
                -o {output} \
                {params.bam}
        """

rule split_BAM:
    input:
        indexlist="data/splitBAM/{seqrun_id}.indexlist.txt",
        bam="data/splitBAM/{seqrun_id}.nsorted.bam"
    output:
        touch("data/splitBAM/{seqrun_id}.done")
    message: "Split BAM based on perfect index matching: {wildcards.seqrun_id}"
    params:
        dir = "data/splitBAM/{seqrun_id}",
        splitBAM = "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/scripts/splitBAM.py"
    log:
        "data/splitBAM/{seqrun_id}_splitBAM.log"
    shell:
        """
        mkdir -p {params.dir}
        python {params.splitBAM} \
                -i {input.bam} \
                -l {input.indexlist} \
                -o {params.dir} > {log}
        """
