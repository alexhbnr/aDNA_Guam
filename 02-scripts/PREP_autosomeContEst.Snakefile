################################################################################
# Project: Ancient human remains from Guam
# Part: Preparation
# Step: Merge autosomal sequencing data per library for contamination estimates
#
# I will merge all sequencing runs of a single library for the contamination
# estimates based on autosomal data.
#
# Alex Huebner, 17/11/19
################################################################################

from glob import glob
import os.path
import re

import pysam

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam"

#### SAMPLES ###################################################################
extract_libid = re.compile("(SP[0-9]+-[A-Z][0-9]+)-[0-9]+_[A-Z]+[0-9]+_[0-9]+_lane[0-9].+")
# Specify males based on ratio of X to autosome sequencing data
MALES = ['SP4210']
# List all seqruns of the male samples with autosome capture data
SEQRUNS = [os.path.basename(fn).replace(".bam", "")
           for sm in MALES
           for fn in glob(f"analysis/{sm}/{sm}*.bam")
           if 'uniq' in fn and len(pysam.AlignmentFile(fn).header['SQ']) > 1]
# Identify all runs belonging to same sample and library
LIBS = {}
for seqrun in SEQRUNS:
    libid = extract_libid.search(seqrun).group(1)
    if "deam" in seqrun:
        flt = "deam"
    else:
        flt = "all"
    if libid + "_" + flt not in LIBS:
        LIBS[libid + "_" + flt] = [seqrun]
    else:
        LIBS[libid + "_" + flt].append(seqrun)
################################################################################


rule all:
    input: expand("analysis/tmp/bam_per_lib/{libid}.bam", libid=LIBS.keys())


rule merge_bams:
    output: "analysis/tmp/bam_per_lib/{libid}.bam"
    message: "Merge the files of sample {wildcards.libid}"
    params: 
        bams = lambda wildcards: " ".join([f"analysis/{libid[:6]}/{libid}.bam" for libid in LIBS[wildcards.libid]])
    shell:
        """
        samtools cat {params.bams} | \
        samtools sort -o {output} -
        samtools index {output}
        """
