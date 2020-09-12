################################################################################
# Project: Ancient human remains from Guam
# Part: Quality
# Step: Determine sex using shotgun data
#
# Alex Huebner, 12/09/20
################################################################################

from glob import glob
import os
import re

import pandas as pd
import pysam

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/"

#### Samples ###################################################################
SAMPLES = ['SP4210', 'SP4211']
# Identify shotgun runs
ARRAYS = pd.read_csv("documentation/capture_arrays.tsv", sep="\t", index_col=[0])
SHOTGUN_RUNS = ARRAYS.loc[ARRAYS['array'] == "shotgun"].index.tolist()
# Identify BAM files 
extract_run = re.compile(r'SP[0-9]+-[A-Z][0-9]+-([0-9]+_[SNDM]+[0-9]+_[0-9]+_lane[0-9])')
LIBS = [os.path.basename(fn)[:-4]
        for sample in SAMPLES
        for fn in glob(f"analysis/{sample}/*.bam")
        if "L35MQ25" not in fn]
SHOTGUN_LIBS = {lib: lib[:6] for lib in LIBS if extract_run.search(lib).group(1) in SHOTGUN_RUNS}
SHOTGUN_LIBS_PERSAMPLE = {}
for lib, sample in SHOTGUN_LIBS.items():
    if sample not in SHOTGUN_LIBS_PERSAMPLE:
        SHOTGUN_LIBS_PERSAMPLE[sample] = [lib]
    else:
        SHOTGUN_LIBS_PERSAMPLE[sample].append(lib)
FLTS = {'all': 'uniq.L35MQ25',
        'deam': 'uniq.L35MQ25.deam'}
################################################################################

rule all:
    input: 
        "results/QUAL_observed_ratio_X_to_autosomes_shotgun.csv"       

rule concatenate_runs:
    output:
        "/mnt/scratch/alexh/tmp/ancGuam/sex_determination/{sample}.{flt}.bam"  
    message: "Concatenate all shotgun runs of sample {wildcards.sample} for filter {wildcards.flt}"
    params:
        runs = lambda wildcards: f"/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/analysis/{wildcards.sample}/{{{','.join(SHOTGUN_LIBS_PERSAMPLE[wildcards.sample])}}}.{FLTS[wildcards.flt]}.bam",
    shell:
        """
        samtools cat {params.runs} | \
        samtools sort -o {output} -
        samtools index {output}
        """

rule pileup_on_sites_captured_by_array:
    input:
        "/mnt/scratch/alexh/tmp/ancGuam/sex_determination/{sample}.{flt}.bam"  
    output:
        "/mnt/scratch/alexh/tmp/ancGuam/sex_determination/{sample}.{flt}.pileup.gz"
    message: "Generate pileup for {wildcards.sample} for filter {wildcards.flt}"
    shell:
        """
        samtools mpileup {input} | gzip > {output}
        """

rule ratio:
    input:
        bam = "/mnt/scratch/alexh/tmp/ancGuam/sex_determination/{sample}.{flt}.bam",
        pileup = "/mnt/scratch/alexh/tmp/ancGuam/sex_determination/{sample}.{flt}.pileup.gz"
    output:
        "analysis/qual/sex/{sample}.{flt}.ratio.txt"
    message: "Determine the ratio of X to autosome coverage for sample {wildcards.sample} for filter {wildcards.flt}"
    run:
        hg19_chrlengths = pd.DataFrame([list(pysam.AlignmentFile(input.bam).header.to_dict()['SQ'][i].values())
                                        for i in range(0, 23)], columns=['chrom', 'length']) \
                            .set_index(['chrom'])

        pileup = pd.read_csv(input.pileup, sep="\t",
                             index_col=[0], usecols=[0, 1, 3],
                             header=None, names=['chrom', 'pos', 'depth'],
                             dtype={'chrom': str})
        depth = pileup[['depth']].groupby(['chrom']).agg(sum)
        depth = depth.loc[[str(i) for i in range(1, 23)] + ['X']]
        depth['autosomes'] = ~depth.index.isin(['X'])
        depth = depth.join(hg19_chrlengths)
        depth['ratio'] = depth['depth'] / depth['length']

        obs_female = depth.groupby(['autosomes']).agg({'ratio': 'mean'}). \
            transpose().reset_index(drop=True)
        obs_female.columns = ['X', 'autosomes']
        obs_female['ratio'] = obs_female['X'] / obs_female['autosomes']
        obs_female.assign(sample=wildcards.sample) \
                  .assign(readType=wildcards.flt)[['sample', 'readType', 'X', 'autosomes', 'ratio']] \
            .to_csv(output[0], sep="\t", index=False)

rule summary:
    input:
        expand("analysis/qual/sex/{sample}.{flt}.ratio.txt", sample=SAMPLES, flt=FLTS.keys())
    output:
        "results/QUAL_observed_ratio_X_to_autosomes_shotgun.csv"
    message: "Summarise the X to autosome ratios"
    run:
        pd.concat([pd.read_csv(fn, sep="\t") for fn in input]) \
            .to_csv(output[0], sep="\t", index=False)
