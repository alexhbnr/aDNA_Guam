################################################################################
# Project: Ancient human remains from Guam
# Part: Data preparation
# Step: Genotyping
#
# We will determine pseudo-haploid genotypes for all sites covered by the
# sequencing data of a sample. We will apply to the following strategy:
#
# We will use an approach used by Mateja Hajdinjak at MPI EVA. We will
# mask terminal C to T substitutions on the terminal 5 sites of a read and
# perform random allele sampling on them.
#
# Alex Huebner, 17/11/19
################################################################################

from glob import glob
import os.path
import re

import pandas as pd

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam"

#### Samples ###################################################################
SAMPLES, = glob_wildcards("analysis/bam/{sample}.nuclear.bam")
################################################################################

#### Genotype info #############################################################
FLT = ['all', 'deam']
CHROMS = [str(i) for i in range(1, 23)] + ['X', 'Y']
################################################################################

wildcard_constraints:
    sample = "SP[0-9]+[A-Z]*\.(all|deam)",
    array = "(1240K|archaicAdmix)",
    trim = "(trim|nontrim)",
    chr = "[0-9XY]+"

rule all:
    input:
        expand("analysis/genotypes/{sample}.vcf.gz", sample=SAMPLES)

rule maskTerminalDeam:
    # Replace Ts in the terminal 5 bases on each end of the read with Ns
    output:
        bam = temp("analysis/genotypes/{sample}.sorted.bam"),
        bai = temp("analysis/genotypes/{sample}.sorted.bam.bai")
    message: "Replace the 5 terminal bases on each end of sample {wildcards.sample} with Ns"
    params:
        maskTerminalDeam = "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/scripts/maskTerminalDeam.py",
        bam = "analysis/bam/{sample}.nuclear.bam",
        reffa = "/mnt/solexa/Genomes/hg19_evan/whole_genome.fa"
    shell:
        """
        {params.maskTerminalDeam} \
            --input {params.bam} \
            --left 5 \
            --right 5 \
            --output /dev/stdout | \
        samtools calmd -uQ - {params.reffa} | \
        samtools sort -o {output.bam} -
        samtools index {output.bam}
        """

rule bamsample_maskedTerminalDeam:
    # Random-allele sampling
    input:
        bam = "analysis/genotypes/{sample}.sorted.bam",
        bai = "analysis/genotypes/{sample}.sorted.bam.bai"
    output:
        vcf = "analysis/genotypes/{sample}.vcf.gz",
        tbi = "analysis/genotypes/{sample}.vcf.gz.tbi"
    message: "Random allele sampling from the reads of sample {wildcards.sample}"
    params:
        bamcaller = "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/scripts/bam-caller.py",
        outputprefix = "analysis/genotypes/{sample}"
    shell:
        """
        if [ -f {params.outputprefix}.vcf ]; then
            rm {params.outputprefix}.vcf
        fi
        {params.bamcaller} \
                --bam {input.bam} \
                --strategy random \
                --seed 0 \
                --mincov 1 \
                --minbq 30 \
                --minmq 25 \
                --sample-name {wildcards.sample} \
                --output {params.outputprefix}
        bgzip {params.outputprefix}.vcf
        tabix {params.outputprefix}.vcf.gz
        """
################################################################################
