################################################################################
# Project: Ancient human remains from Guam
# Part: Quality
# Step: Determine sex
#
# Alex Huebner, 17/11/19
################################################################################

from functools import reduce
from glob import glob
import os.path
import re
import sys

import pysam
import numpy as np
import pandas as pd

from snakemake.utils import R

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/"

#### Samples ###################################################################
SAMPLES = ['SP4210', 'SP4211']
# BAMS
LIBS = {os.path.basename(fn)[:-4]: sample
        for sample in SAMPLES
        for fn in glob(f"analysis/{sample}/*.bam") 
        if "L35MQ25" not in fn}
# CAPTURE ARRAY
ARRAYS = pd.read_csv("documentation/capture_arrays.tsv", sep="\t", index_col=[0])
LIBS = {lib: sample for lib, sample in LIBS.items() if ARRAYS.loc[lib.split("-")[2], 'array'] in ['390K', '390Ksupp']}
# CAPTURE ARRAY SITES
ARRAY_BEDS = {'390K': "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/documentation/390K.bed.gz",
              '390Ksupp': "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/documentation/840K.bed.gz"}
# Filters
FLTS = {'all': 'uniq.L35MQ25',
        'deam': 'uniq.L35MQ25.deam'}
################################################################################

rule all:
    input:
        "analysis/qual/sex/sex.csv"

rule pileup_on_sites_captured_by_array:
    output:
        temp("/mnt/scratch/alexh/tmp/ancFlores/sex_determination/{lib}.{flt}.pileup.gz")
    message: "Generate pileup on subset of the sites captured by the array for {wildcards.lib}"
    params:
        bam = lambda wildcards: f"analysis/{LIBS[wildcards.lib]}/{wildcards.lib}.{FLTS[wildcards.flt]}.bam",
        bed = lambda wildcards: ARRAY_BEDS[ARRAYS.at[wildcards.lib.split("-")[2], 'array']]
    shell:
        """
        samtools mpileup \
            -l {params.bed} {params.bam} | gzip > {output}
        """

rule cummulative_coverage:
    input:
        "/mnt/scratch/alexh/tmp/ancFlores/sex_determination/{lib}.{flt}.pileup.gz"
    output:
        "analysis/qual/sex/{lib}.{flt}.bp_per_chr.txt"
    message: "Calculate the cummulative coverage per chromosome for library {wildcards.lib} filtered for {wildcards.flt} reads"
    run:
        pileup = pd.read_csv(input[0],
                             sep="\t", index_col=[0], usecols=[0, 1, 3],
                             header=None, names=['chrom', 'pos', 'depth'],
                             dtype={'chrom': str})
        pileup[['depth']].groupby(['chrom']).agg(sum). \
            to_csv(output[0], sep="\t")

rule sex_determination:
    input:
        expand("analysis/qual/sex/{lib}.{flt}.bp_per_chr.txt", lib=LIBS.keys(), flt=FLTS.keys())
    output:
        "analysis/qual/sex/sex.csv"
    message: "Summarise per SNP coverages and determine sex per sample"
    run:
        def calculate_ratio(df, sample, array, flt):
            if 'X' not in df.index:
                df.loc['X'] = [0]
            if 'Y' in df.index:
                df = df.drop(['Y'])
            df['autosomes'] = ~df.index.isin(['X', 'Y'])
            obs_female = df.groupby(['autosomes']).agg(['mean']). \
                transpose().reset_index().drop(['level_1'], axis=1)
            obs_female.columns = ['array', 'X', 'autosomes']
            obs_female['observed'] = obs_female['X'] / obs_female['autosomes']
            obs_female['array'] = array
            obs_female['sample'] = sample
            obs_female['flt'] = flt
            return obs_female

        # Number of sites per capture-array and chromosome
        arrays = [pd.read_csv(ARRAY_BEDS[array], sep="\t", index_col=[0],
                              header=None, usecols=[0, 1],
                              names=['chrom', array],
                              dtype={'chrom': str}). \
                  groupby(level=0).count()
                  for array in ARRAY_BEDS]
        capture_sites = reduce(lambda left, right: left.join(right), arrays)
        # Expected ratio for females
        capture_sites['autosomes'] = ~capture_sites.index.isin(['X', 'Y'])
        capture_sites = capture_sites.drop(['Y'])
        expected_female = capture_sites.groupby(['autosomes']).agg(['mean']). \
            transpose().reset_index().drop(['level_1'], axis=1)
        expected_female.columns = ['array', 'X', 'autosomes']
        expected_female['ratio'] = expected_female['X'] / expected_female['autosomes']
        expected_female = expected_female.set_index(['array'])

        # Number of alleles observed per library
        libs = [pd.read_csv(fn, sep="\t", index_col=[0],
                            dtype={'chrom': str})
                for fn in glob("analysis/qual/sex/*.bp_per_chr.txt")]
        extract_info = re.compile(r'(SP[0-9]+)-[A-Z][0-9]+-([0-9]+_[A-Z]+[0-9]+_[0-9]+_lane[0-9])\.([a-z]+).bp_per_chr.txt')
        lib_attr = [extract_info.search(fn).groups()
                    for fn in glob("analysis/qual/sex/*.bp_per_chr.txt")]
        lib_ratios = pd.concat([calculate_ratio(lib, attrs[0], ARRAYS.loc[attrs[1], 'array'], attrs[2])
                                for lib, attrs in zip(libs, lib_attr)])
        lib_ratios_summary = lib_ratios[['sample', 'flt', 'array', 'observed']]. \
                                groupby(['sample', 'flt', 'array']).agg('mean'). \
                             reset_index()
        lib_ratios_summary['expected'] = expected_female.loc[lib_ratios_summary['array'], 'ratio'].values
        lib_ratios_summary['ratio'] = lib_ratios_summary['observed'] / lib_ratios_summary['expected']
        lib_ratios_summary['sex'] = ["male" if r < 0.6 else "female" for r in lib_ratios_summary['ratio']]
        lib_ratios_summary.to_csv(output[0], sep="\t", index=False)

################################################################################
