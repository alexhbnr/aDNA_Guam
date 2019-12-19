################################################################################
# Project: Ancient human remains from Guam
# Part: Preparation
# Step: Merge all autosomal capture and all MT capture sequencing data per
#       sample
#
# Alex Huebner, 17/11/19
################################################################################

from glob import glob
from os.path import basename
import re

import pysam

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam"

SAMPLES = ['SP4210', 'SP4211']

FILEEXTS = {'all': 'uniq.L35MQ25',
            'deam': 'uniq.L35MQ25.deam'}

def list_autosome_bamfiles(sample, flt):
    '''Generate a list of all BAM files that contain autosomal capture data for each
       filtering step
    '''
    extract_libid = re.compile("SP[0-9]+-([A-Z][0-9]+)-[0-9]+_[A-Z]+[0-9]+_[0-9]+_lane[0-9]")
    if sample == "SP4163HQ":
        sample = "SP4163"
        seqruns = [seqrun
                   for seqrun in glob(f"analysis/{sample}/{sample}-*.{flt}.bam")
                   if (len(pysam.AlignmentFile(seqrun).header['SQ']) > 1 and
                       extract_libid.search(basename(seqrun).replace(f".{flt}.bam", "")).group(1)) in SP4163HQLIBS]
    else:
        seqruns = [seqrun
                   for seqrun in glob(f"analysis/{sample}/{sample}-*.{flt}.bam")
                   if len(pysam.AlignmentFile(seqrun).header['SQ']) > 1]
    return " ".join(seqruns)


def list_MT_bamfiles(sample, flt):
    '''Generate a list of all BAM files that contain MT capture data for each
       filtering step
    '''
    seqruns = [seqrun
               for seqrun in glob(f"analysis/{sample}/{sample}-*.{flt}.bam")
               if len(pysam.AlignmentFile(seqrun).header['SQ']) == 1]
    return " ".join(seqruns)


def extract_libids_MT(sample, flt):
    '''Generate a list of all lib ids that were captured for MT
    '''
    extract_libid = re.compile("SP[0-9]+-([A-Z][0-9]+)-[0-9]+_[A-Z]+[0-9]+_[0-9]+_lane[0-9]")
    libids = [extract_libid.search(basename(seqrun).replace(f".{flt}.bam", "")).group(1)
              for seqrun in glob(f"analysis/{sample}/{sample}-*.{flt}.bam")
              if len(pysam.AlignmentFile(seqrun).header['SQ']) == 1]
    libids = list(set(libids))
    return [{"ID": libid, "SM": sample} for libid in libids]


def extract_libids_autosomes(sample, flt):
    '''Generate a list of all lib ids that were captured using the 1240K array
    '''
    extract_libid = re.compile("SP[0-9]+-([A-Z][0-9]+)-[0-9]+_[A-Z]+[0-9]+_[0-9]+_lane[0-9]")
    if sample == "SP4163HQ":
        sample = "SP4163"
        libids = [extract_libid.search(basename(seqrun).replace(f".{flt}.bam", "")).group(1)
                  for seqrun in glob(f"analysis/{sample}/{sample}-*.{flt}.bam")
                  if len(pysam.AlignmentFile(seqrun).header['SQ']) > 1]
        libids = [lid for lid in libids if lid in SP4163HQLIBS]
    else:
        libids = [extract_libid.search(basename(seqrun).replace(f".{flt}.bam", "")).group(1)
                  for seqrun in glob(f"analysis/{sample}/{sample}-*.{flt}.bam")
                  if len(pysam.AlignmentFile(seqrun).header['SQ']) > 1]
    libids = list(set(libids))
    return [{"ID": libid, "SM": sample} for libid in libids]


rule all:
    input:  expand("analysis/bam/{sample}.{flt}.nuclear.bam", sample=SAMPLES, flt=['all', 'deam']),
            expand("analysis/bam/{sample}.{flt}.MT.bam", sample=[sample for sample in SAMPLES if sample != "SP4163HQ"], flt=['all', 'deam'])

rule prepare_header_MT:
    output: 
        temp("analysis/bam/{sample}.{flt}.MT.header")
    message: "Prepare BAM header for sample {wildcards.sample} for MT capture data"
    params:
        rgs = lambda wildcards: extract_libids_MT(wildcards.sample, FILEEXTS[wildcards.flt])
    run:
        header = {'HD': {'VN': '1.4'},
                  'SQ': [{'LN': 16569, 'SN': 'MT'}],
                 }
        header['RG'] = params.rgs
        with pysam.AlignmentFile(output[0], "wh", header=header) as outfile:
            outfile.close()


rule prepare_header_autosomes:
    output:
        temp("analysis/bam/{sample}.{flt}.nuclear.header")
    message: "Prepare BAM header for sample {wildcards.sample} for autosome capture data"
    params:
        rgs = lambda wildcards: extract_libids_autosomes(wildcards.sample, FILEEXTS[wildcards.flt])
    run:
        header = {'HD': {'VN': '1.4'},
                  'SQ': [{'LN': 249250621, 'SN': '1'},
                         {'LN': 243199373, 'SN': '2'},
                         {'LN': 198022430, 'SN': '3'},
                         {'LN': 191154276, 'SN': '4'},
                         {'LN': 180915260, 'SN': '5'},
                         {'LN': 171115067, 'SN': '6'},
                         {'LN': 159138663, 'SN': '7'},
                         {'LN': 146364022, 'SN': '8'},
                         {'LN': 141213431, 'SN': '9'},
                         {'LN': 135534747, 'SN': '10'},
                         {'LN': 135006516, 'SN': '11'},
                         {'LN': 133851895, 'SN': '12'},
                         {'LN': 115169878, 'SN': '13'},
                         {'LN': 107349540, 'SN': '14'},
                         {'LN': 102531392, 'SN': '15'},
                         {'LN': 90354753, 'SN': '16'},
                         {'LN': 81195210, 'SN': '17'},
                         {'LN': 78077248, 'SN': '18'},
                         {'LN': 59128983, 'SN': '19'},
                         {'LN': 63025520, 'SN': '20'},
                         {'LN': 48129895, 'SN': '21'},
                         {'LN': 51304566, 'SN': '22'},
                         {'LN': 155270560, 'SN': 'X'},
                         {'LN': 59373566, 'SN': 'Y'},
                         {'LN': 17569, 'SN': 'MT'},
                         {'LN': 4262, 'SN': 'GL000207.1'},
                         {'LN': 15008, 'SN': 'GL000226.1'},
                         {'LN': 19913, 'SN': 'GL000229.1'},
                         {'LN': 27386, 'SN': 'GL000231.1'},
                         {'LN': 27682, 'SN': 'GL000210.1'},
                         {'LN': 33824, 'SN': 'GL000239.1'},
                         {'LN': 34474, 'SN': 'GL000235.1'},
                         {'LN': 36148, 'SN': 'GL000201.1'},
                         {'LN': 36422, 'SN': 'GL000247.1'},
                         {'LN': 36651, 'SN': 'GL000245.1'},
                         {'LN': 37175, 'SN': 'GL000197.1'},
                         {'LN': 37498, 'SN': 'GL000203.1'},
                         {'LN': 38154, 'SN': 'GL000246.1'},
                         {'LN': 38502, 'SN': 'GL000249.1'},
                         {'LN': 38914, 'SN': 'GL000196.1'},
                         {'LN': 39786, 'SN': 'GL000248.1'},
                         {'LN': 39929, 'SN': 'GL000244.1'},
                         {'LN': 39939, 'SN': 'GL000238.1'},
                         {'LN': 40103, 'SN': 'GL000202.1'},
                         {'LN': 40531, 'SN': 'GL000234.1'},
                         {'LN': 40652, 'SN': 'GL000232.1'},
                         {'LN': 41001, 'SN': 'GL000206.1'},
                         {'LN': 41933, 'SN': 'GL000240.1'},
                         {'LN': 41934, 'SN': 'GL000236.1'},
                         {'LN': 42152, 'SN': 'GL000241.1'},
                         {'LN': 43341, 'SN': 'GL000243.1'},
                         {'LN': 43523, 'SN': 'GL000242.1'},
                         {'LN': 43691, 'SN': 'GL000230.1'},
                         {'LN': 45867, 'SN': 'GL000237.1'},
                         {'LN': 45941, 'SN': 'GL000233.1'},
                         {'LN': 81310, 'SN': 'GL000204.1'},
                         {'LN': 90085, 'SN': 'GL000198.1'},
                         {'LN': 92689, 'SN': 'GL000208.1'},
                         {'LN': 106433, 'SN': 'GL000191.1'},
                         {'LN': 128374, 'SN': 'GL000227.1'},
                         {'LN': 129120, 'SN': 'GL000228.1'},
                         {'LN': 137718, 'SN': 'GL000214.1'},
                         {'LN': 155397, 'SN': 'GL000221.1'},
                         {'LN': 159169, 'SN': 'GL000209.1'},
                         {'LN': 161147, 'SN': 'GL000218.1'},
                         {'LN': 161802, 'SN': 'GL000220.1'},
                         {'LN': 164239, 'SN': 'GL000213.1'},
                         {'LN': 166566, 'SN': 'GL000211.1'},
                         {'LN': 169874, 'SN': 'GL000199.1'},
                         {'LN': 172149, 'SN': 'GL000217.1'},
                         {'LN': 172294, 'SN': 'GL000216.1'},
                         {'LN': 172545, 'SN': 'GL000215.1'},
                         {'LN': 174588, 'SN': 'GL000205.1'},
                         {'LN': 179198, 'SN': 'GL000219.1'},
                         {'LN': 179693, 'SN': 'GL000224.1'},
                         {'LN': 180455, 'SN': 'GL000223.1'},
                         {'LN': 182896, 'SN': 'GL000195.1'},
                         {'LN': 186858, 'SN': 'GL000212.1'},
                         {'LN': 186861, 'SN': 'GL000222.1'},
                         {'LN': 187035, 'SN': 'GL000200.1'},
                         {'LN': 189789, 'SN': 'GL000193.1'},
                         {'LN': 191469, 'SN': 'GL000194.1'},
                         {'LN': 211173, 'SN': 'GL000225.1'},
                         {'LN': 547496, 'SN': 'GL000192.1'},
                         {'LN': 171823, 'SN': 'NC_007605'},
                         {'LN': 35477943, 'SN': 'hs37d5'},
                         {'LN': 6386, 'SN': 'phiX'}]
                 }
        header['RG'] = params.rgs
        with pysam.AlignmentFile(output[0], "wh", header=header) as outfile:
            outfile.close()

rule merge_BAM_autosomes:
    input:
        "analysis/bam/{sample}.{flt}.nuclear.header"
    output:
        "analysis/bam/{sample}.{flt}.nuclear.bam"
    message: "Merge all sequencing files of {wildcards.sample} filtered for {wildcards.flt} reads of the autosome capture array"
    params:
        bams = lambda wildcards: list_autosome_bamfiles(wildcards.sample, FILEEXTS[wildcards.flt]),
        reffa = "/mnt/solexa/Genomes/hg19_evan/whole_genome.fa"
    shell:
        """
        samtools cat -h {input} {params.bams} | \
        samtools calmd -Q - {params.reffa} | \
        samtools sort -o {output} -
        samtools index {output}
        """


rule merge_BAM_MT:
    input:
        "analysis/bam/{sample}.{flt}.MT.header"
    output:
        "analysis/bam/{sample}.{flt}.MT.bam"
    message: "Merge all sequencing files of {wildcards.sample} filtered for {wildcards.flt} reads captured for MT array"
    params:
        bams = lambda wildcards: list_MT_bamfiles(wildcards.sample, FILEEXTS[wildcards.flt])
    shell:
        """
        samtools cat -h {input} {params.bams} | \
        samtools sort -o {output} -
        samtools index {output}
        """
