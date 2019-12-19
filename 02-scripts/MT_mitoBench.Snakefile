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

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/"

rule all:
    input: "analysis/mtDNA/mitoBench_perseqfile.config"

rule link_plus_samplelist:
    output:
        "analysis/mtDNA/samplelist.txt"
    message: "Identify mtDNA sequencing data, link with shorter name and add to samplelist"
    params:
        bamdir = '/mnt/scratch/alexh/tmp/ancGuam/mitoBench_persample/bamfiles',
    run:
        cwd = getcwd()
        # Prepare output directory
        makedirs(params.bamdir, exist_ok=True)
        # Prepare sample list
        samplefile = open(output[0], "wt")

        # Identify all Flores samples
        for fn in glob(f"{cwd}/analysis/bam/*.MT.bam"):
            sample = basename(fn).replace(".MT.bam", "")
            samplefile.write(f"{sample}\n") 
            symlink(fn, f"{params.bamdir}/{sample}.bam")

        samplefile.close()

rule write_configfile:
    input:
        "analysis/mtDNA/samplelist.txt"
    output:
        "analysis/mtDNA/mitoBench_perseqfile.config"
    message: "Generate the config file to run the mitoBench ancient mtDNA pipeline on single sequencing data files"
    params:
        bamdir = '/mnt/scratch/alexh/tmp/ancGuam/mitoBench_persample/bamfiles',
        projdir = 'analysis/mtDNA',
        tmpdir = '/mnt/scratch/alexh/tmp/ancGuam/mitoBench_persample'
    run:
        cwd = getcwd()
        with open(output[0], "wt") as outfile:
            outfile.write('{\n')
            outfile.write(f'    "samplelist" : "{cwd}/{input[0]}",\n')
            outfile.write(f'    "bamdir" : "{params.bamdir}",\n')
            outfile.write(f'    "projdir" : "{cwd}/{params.projdir}",\n')
            outfile.write(f'    "tmpdir" : "{params.tmpdir}",\n')
            outfile.write('    "bamsuffix" : "bam",\n')
            outfile.write('    "sampleIDconstraint" : "SP[0-9]+.[a-z]+"\n')
            outfile.write('}\n')
