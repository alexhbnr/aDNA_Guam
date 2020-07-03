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
              '390Ksupp': "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/documentation/840K.bed.gz",
              '1240K': "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/documentation/390Kplus840K.bed.gz"}
# Filters
FLTS = {'all': 'uniq.L35MQ25',
        'deam': 'uniq.L35MQ25.deam'}
################################################################################

rule all:
    input:
        f"{workflow.basedir}/../05-results/QUAL_observed_ratio_X_to_autosomes.csv",
        f"{workflow.basedir}/../05-results/QUAL_expected_ratio_X_to_autosomes.csv"

rule pileup_on_sites_captured_by_array:
    output:
        temp("/mnt/scratch/alexh/tmp/ancFlores/sex_determination/{lib}.{flt}.pileup.gz")
    message: "Generate pileup on subset of the sites captured by the array for {wildcards.lib}"
    params:
        bam = lambda wildcards: f"analysis/{LIBS[wildcards.lib]}/{wildcards.lib}.{FLTS[wildcards.flt]}.bam",
        bed = lambda wildcards: ARRAY_BEDS['1240K']
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
        obs = f"{workflow.basedir}/../05-results/QUAL_observed_ratio_X_to_autosomes.csv",
        exp = f"{workflow.basedir}/../05-results/QUAL_expected_ratio_X_to_autosomes.csv"
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
        capture_sites['1240K'] = capture_sites['390K'] + capture_sites['390Ksupp']
        # Expected ratio for females
        capture_sites['autosomes'] = ~capture_sites.index.isin(['X', 'Y'])
        capture_sites = capture_sites.drop(['Y'])
        expected_female = pd.DataFrame.from_dict({'autosome': [capture_sites.loc[capture_sites.index != 'X', '1240K'].sum()],
                                                  'X': [capture_sites.loc[capture_sites.index == 'X', '1240K'][0]]})
        expected_female['female'] = expected_female['X'] / (expected_female['X'] + expected_female['autosome'])
        expected_female['male'] = expected_female['X'] / (expected_female['X'] + 2 * expected_female['autosome'])

        # Number of alleles observed per library
        extract_info = re.compile(r'(SP[0-9]+)-([A-Z][0-9]+)-([0-9]+_[A-Z]+[0-9]+_[0-9]+_lane[0-9])\.([a-z]+).bp_per_chr.txt')
        lib_attr = pd.DataFrame([extract_info.search(f"analysis/qual/sex/{lib}.{flt}.bp_per_chr.txt").groups()
                                 for lib in LIBS.keys()
                                 for flt in FLTS.keys()],
                                 columns=['sample', 'library', 'run', 'flt'])
        lib_attr['prefix'] = [f'{lib}.{flt}' for lib in LIBS.keys() for flt in FLTS.keys()]
        lib_attr['array'] = [ARRAYS.loc[run, 'array'] for run in lib_attr['run']]

        libs = [pd.read_csv(fn, sep="\t", index_col=[0],
                            dtype={'chrom': str})
                for fn in glob("analysis/qual/sex/*.bp_per_chr.txt")]
        # Infer the array state of each run
        def calculate_ratio(runs):
            index = [str(i) for i in range(1, 23)] + ['X', 'Y']
            # Calculate cumulative coverage across all libraries filter by the same scheme
            cumcov = pd.read_csv(f'analysis/qual/sex/{runs[0]}.bp_per_chr.txt',
                                 sep="\t", index_col=[0]) \
                .reindex(index) \
                .drop(['Y'])
            for i, run in enumerate(runs[1:]):
                cov_run = pd.read_csv(f'analysis/qual/sex/{run}.bp_per_chr.txt',
                                      sep="\t", index_col=[0]) \
                    .reindex(index) \
                    .drop(['Y'])
                cov_run = cov_run.fillna(value=0)
                cumcov += cov_run

            obs_female = pd.DataFrame.from_dict({'autosome': [cumcov.loc[cumcov.index != 'X', 'depth'].sum()],
                                                 'X': [cumcov.loc[cumcov.index == 'X', 'depth'][0]]})
            obs_female['female'] = obs_female['X'] / (obs_female['X'] + obs_female['autosome'])
            return obs_female

        x_auto_ratio = lib_attr.groupby(['sample', 'flt']).apply(lambda x: calculate_ratio(x['prefix'].tolist())) \
            .reset_index().drop(['level_2'], axis=1) \
            .to_csv(output.obs, sep="\t", index=False)
        expected_female.to_csv(output.exp, sep="\t", index=False)

################################################################################
