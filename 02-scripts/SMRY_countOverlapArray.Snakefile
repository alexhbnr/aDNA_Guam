################################################################################
# Project: Ancient human remains from Guam
# Part: Summary
# Step: Count the number of SNPs of each nuclear array covered by the
#       sequencing data
#
# Alex Huebner, 18/11/19
################################################################################

from glob import glob
import os.path
import re

import pandas as pd

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam"

#### Samples ###################################################################
SAMPLES = ['SP4210', 'SP4211']
# Capture array info - skip MT
CAPTURE_ARRAY = pd.read_csv(f"{workflow.basedir}/../01-documentation/capture_arrays.tsv", sep="\t", index_col=[0])
extract_seqrun = re.compile("SP[0-9]+-[A-Z][0-9]+-([0-9]+_[A-Z]+[0-9]+_[0-9]+_lane[0-9])\..+")
BAMS = {os.path.basename(fn)[:-4]: sample
        for sample in SAMPLES
        for fn in glob(f"analysis/{sample}/*.bam")
        if ("L35MQ25" in fn and
             CAPTURE_ARRAY.loc[extract_seqrun.search(os.path.basename(fn)).group(1), 'array'] != "MT")}
################################################################################

#### Arrays ####################################################################
ARRAY = {'390K': "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/documentation/390K.bed.gz",
         '390Ksupp': "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/documentation/840K.bed.gz",
         'archaicAdmix': "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/documentation/Archaic.align.noN.sorted.bed.gz"}
################################################################################

wildcard_constraints:
    lib = "SP[0-9]+-[A-Z][0-9]+-[0-9]+_[A-Z]+[0-9]+_[0-9]+_lane[0-9]\.uniq.(L35MQ25|L35MQ25\.deam)",
    array = "(390K|390Ksupp|archaicAdmix)",


rule all:
    input: 
        expand("analysis/logs/nSNPs/{lib}.{array}.nsnps.txt", lib=BAMS.keys(), array=ARRAY.keys())

rule count_snps:
    output:
        "analysis/logs/nSNPs/{lib}.{array}.nsnps.txt"
    message: "Count number of SNPs covered by library {wildcards.lib} for array {wildcards.array}"
    params: 
        bam = lambda wildcards: f"analysis/{BAMS[wildcards.lib]}/{wildcards.lib}.bam",
        tmpbam = "/mnt/scratch/alexh/tmp/{lib}.bam",
        array = lambda wildcards: ARRAY[wildcards.array]
    shell:
        """
        ln -s ${{PWD}}/{params.bam} {params.tmpbam}
        samtools index {params.tmpbam}
        nsnps=$(zcat {params.array} | \
                samtools mpileup \
                    -q 30 \
                    -Q 25 \
                    -l - \
                    {params.tmpbam} | wc -l)
        echo -e "{wildcards.lib}\t{wildcards.array}\t${{nsnps}}" > {output}
        rm {params.tmpbam}*
        """
