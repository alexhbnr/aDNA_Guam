################################################################################
# Project: Ancient human remains from Guam
# Part: Quality
# Step: Prepare files for mitoBench ancient mtDNA pipeline
#
# In order to estimate the contamination from each library, I will run the
# mitoBench ancient mtDNA pipeline on each individual sequencing file. For
# per-sample mtDNA analysis will later be run separately after merging all
# mtDNA-capture data of a sample.
#
# The mitoBench pipeline will be manually executed to avoid having nested
# Snakemake instances.
#
# Alex Huebner, 17/11/19
################################################################################

from glob import glob
from os import getcwd, makedirs
from os.path import basename
import re
import subprocess

import pysam
import pandas as pd

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/"


rule all:
    input: "analysis/qual/mtDNA_contamination/mitoBench_perseqfile.config"

rule link_plus_samplelist:
    output:
        "analysis/qual/mtDNA_contamination/samplelist.txt"
    message: "Identify mtDNA sequencing data, link with shorter name and add to samplelist"
    params:
        bamdir = '/mnt/scratch/alexh/tmp/ancGuam/mitoBench/bamfiles',
    run:
        # Extract library id for shorter file names
        extract_libid = re.compile("SP[0-9]+-([A-Z][0-9]+)-[0-9]+_[A-Z]+[0-9]+_[0-9]+_lane[0-9].bam")
        # Prepare output directory
        makedirs(params.bamdir, exist_ok=True)
        # Prepare sample list
        samplefile = open(output[0], "wt")

        # Identify all Flores samples
        samples = [basename(fld) for fld in glob("analysis/*")
                   if basename(fld).startswith("SP")]
        for sample in samples:  # Identify all sequencing data files of a sample
            seqruns = [seqrun for seqrun in glob(f'analysis/{sample}/*.bam') if 'uniq' not in seqrun]
            for seqrun in seqruns:  # Identify if sequencing data is aligned against MT genome
                ref_contigs = pysam.AlignmentFile(seqrun).header['SQ']
                if (len(ref_contigs) == 1 and
                    ref_contigs[0]['SN'] == 'gi|251831106|ref|NC_012920.1|'):
                    cwd = getcwd()
                    libid = extract_libid.search(basename(seqrun)).group(1)
                    seqrun_prefix = seqrun.replace(".bam", "")
                    # All reads
                    samplefile.write(f'{sample}-{libid}_all\tPAIRED\n')
                    subprocess.run(f'samtools view -H {seqrun_prefix}.uniq.L35MQ25.bam | ' +
                                   f'sed "s/gi|251831106|ref|NC_012920.1|/MT/g" | ' +
                                   f'samtools reheader - {seqrun_prefix}.uniq.L35MQ25.bam > ' +
                                   f'{params.bamdir}/{sample}-{libid}_all.bam', shell=True)
                    # Deaminated only reads
                    samplefile.write(f'{sample}-{libid}_deam\n')
                    subprocess.run(f'samtools view -H {seqrun_prefix}.uniq.L35MQ25.deam.bam | ' +
                                   f'sed "s/gi|251831106|ref|NC_012920.1|/MT/g" | ' +
                                   f'samtools reheader - {seqrun_prefix}.uniq.L35MQ25.deam.bam > ' +
                                   f'{params.bamdir}/{sample}-{libid}_deam.bam', shell=True)

        samplefile.close()

rule write_configfile:
    input:
        "analysis/qual/mtDNA_contamination/samplelist.txt"
    output:
        "analysis/qual/mtDNA_contamination/mitoBench_perseqfile.config"
    message: "Generate the config file to run the mitoBench ancient mtDNA pipeline on single sequencing data files"
    params:
        json_template = "/home/alexander_huebner/github//mitoBench-ancientMT/mitoBench_pipeline-config.json",
        bamdir = '/mnt/scratch/alexh/tmp/ancGuam/mitoBench/bamfiles',
        projdir = '/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/analysis/qual/mtDNA_contamination',
        tmpdir = '/mnt/scratch/alexh/tmp/ancGuam/mitoBench'
    run:
        cwd = getcwd()

        with open(params.json_template) as jsonfile:
            config = json.load(jsonfile)

        config['samplelist'] = f"{cwd}/{input[0]}"
        config['seqdatadir'] = params.bamdir
        config['seqdatatype'] = "bam"
        config['seqdatasuffix'] = "bam"
        config['sampleIDconstraint'] = "SP[0-9]+-[A-Z][0-9]+_[a-z]+"
        config['projdir'] = params.projdir
        config['tmpdir'] = params.tmpdir

        with open(output[0], "wt") as outfile:
            json.dump(config, outfile)

################################################################################

#### Summarise #################################################################

rule summary:
    output:
        "results/QUAL_mitoBench_perlibrary.csv"
    message: "Concatenate per-sample summary tables"
    params:
        dir = '/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/analysis/qual/mtDNA_contamination'
    run:
        summary_tables = pd.concat([pd.read_csv(sample + "/summary_table.csv", sep="\t")
                                    for sample in glob("/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/analysis/qual/mtDNA_contamination/SP*")])
        summary_tables['library'] = summary_tables['sample'].str.extract(r'SP421[01]-([A-Z][0-9]+)_[a-z]+')
        summary_tables['readType'] = summary_tables['sample'].str.extract(r'SP421[01]-[A-Z][0-9]+_([a-z]+)') \
                .replace({'all': 'all reads',
                          'deam': 'deaminated only reads'})
        summary_tables['sample'] = summary_tables['sample'].str.extract(r'(SP421[01])-[A-Z][0-9]+_[a-z]+')
        summary_tables[['sample', 'library', 'readType'] + [col for col in summary_tables.columns if col not in ['sample', 'library', 'readType']]] \
            .to_csv(output[0], sep="\t", index=False)
