################################################################################
# Project: Ancient human remains from Guam
# Part: Mitochondrial DNA analysis
# Step: Prepare files for mitoBench ancient mtDNA pipeline of the samples
#
# Alex Huebner, 17/11/19
################################################################################

from glob import glob
from os import getcwd, makedirs, symlink
from os.path import basename
import re
import subprocess

import pysam
import pandas as pd

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/"

#### Identify libraries per sample #############################################
libraries = pd.read_csv("analysis/qual/mtDNA_contamination/summary.csv", sep="\t")
libraries['contamMix'] = libraries['proportion of contamination DNA (contamMix)'].str.extract(r'([0-9]+\.*[0-9]*)% .+').astype(float)
LIBRARIES = libraries.loc[libraries['contamMix'] < 25]
LIBRARIES['seqrun_id'] = ""
LIBRARIES.loc[LIBRARIES['library'].str[0] == "D", 'seqrun_id'] = "160906_M02279_0022_lane1"
LIBRARIES.loc[LIBRARIES['library'].str[0] == "F", 'seqrun_id'] = "170822_M02279_0141_lane1"
LIBRARIES['readType'] = LIBRARIES['readType'].replace({'all reads': 'all',
                                                       'deaminated only reads': 'deam'})
LIBRARIES['filename'] = LIBRARIES['sample'] + "-" + LIBRARIES['library'] + "-" + LIBRARIES['seqrun_id']
FILEPREFIXES = LIBRARIES.groupby(['sample', 'readType'], as_index=False)['filename'].aggregate(lambda x: ",".join(x))
################################################################################


rule all:
    input:
        "analysis/mtDNA/mitoBench_perseqfile.config",
        expand("/mnt/scratch/alexh/tmp/ancGuam/mitoBench_sample/{sample}_{flt}.bam", sample=['SP4210', 'SP4211'], flt=['all', 'deam'])

rule concat_bam:
    output:
        "/mnt/scratch/alexh/tmp/ancGuam/mitoBench_sample/{sample}_{flt}.bam"
    message: "Concatenate BAM files of sample {wildcards.sample} for filter {wildcards.flt}"
    params:
        bam = lambda wildcards: f"/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/analysis/{wildcards.sample}/{{{FILEPREFIXES.loc[(FILEPREFIXES['sample'] == wildcards.sample) & (FILEPREFIXES['readType'] == wildcards.flt), 'filename'].iloc[0]}}}",
        suffix = lambda wildcards: "uniq.L35MQ25.bam" if wildcards.flt == "all" else "uniq.L35MQ25.deam.bam"
    shell:
        """
        samtools cat {params.bam}.{params.suffix} | \
        samtools sort -o {output} -
        """


rule samplelist:
    output:
        "analysis/mtDNA/samplelist.txt"
    message: "Identify mtDNA sequencing data, link with shorter name and add to samplelist"
    params:
        bamdir = '/mnt/scratch/alexh/tmp/ancGuam/mitoBench_sample',
    run:
        # Prepare sample list
        samplefile = open(output[0], "wt")

        samplefile.write('sample\tseqdatatype\n')
        # Identify all Flores samples
        for fn in glob(f"{params.bamdir}/*.bam"):
            sample = basename(fn).split(".")[0]
            samplefile.write(f"{sample}\tSINGLE\n") 

        samplefile.close()

rule write_configfile:
    input:
        "analysis/mtDNA/samplelist.txt"
    output:
        "analysis/mtDNA/mitoBench_perseqfile.config"
    message: "Generate the config file to run the mitoBench ancient mtDNA pipeline on single sequencing data files"
    params:
        json_template = "/home/alexander_huebner/github//mitoBench-ancientMT/mitoBench_pipeline-config.json",
        bamdir = '/mnt/scratch/alexh/tmp/ancGuam/mitoBench_sample/',
        projdir = 'analysis/mtDNA',
        tmpdir = '/mnt/scratch/alexh/tmp/ancGuam/mitoBench_sample'
    run:
        cwd = getcwd()

        with open(params.json_template) as jsonfile:
            config = json.load(jsonfile)

        config['samplelist'] = f"{cwd}/{input[0]}"
        config['seqdatadir'] = params.bamdir
        config['seqdatatype'] = "bam"
        config['seqdatasuffix'] = "bam"
        config['sampleIDconstraint'] = "SP[0-9]+_[a-z]+"
        config['projdir'] = params.projdir
        config['tmpdir'] = params.tmpdir

        with open(output[0], "wt") as outfile:
            json.dump(config, outfile)

################################################################################

#### Summarise #################################################################

rule summary:
    output:
        "results/QUAL_mitoBench_persample.csv"
    message: "Concatenate per-sample summary tables"
    params:
        dir = '/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/analysis/mtDNA/'
    run:
        summary_tables = pd.concat([pd.read_csv(sample + "/summary_table.csv", sep="\t")
                                    for sample in glob(f"{params.dir}/SP*")])
        summary_tables['readType'] = summary_tables['sample'].str.extract(r'SP421[01]_([a-z]+)') \
                .replace({'all': 'all reads',
                          'deam': 'deaminated only reads'})
        summary_tables['sample'] = summary_tables['sample'].str.extract(r'(SP421[01])_[a-z]+')
        summary_tables[['sample', 'readType'] + [col for col in summary_tables.columns if col not in ['sample', 'readType']]] \
            .to_csv(output[0], sep="\t", index=False)
