################################################################################
# Project: Ancient human remains from Guam
# Part: Summary
# Step: Summary of the major information about the nuclear sequencing data
#
# For the publication, I will merge the information on the nuclear sequencing
# data across all libraries of a sample and re-calculate the contamination
# estimates across the merged data.
#
# Alex Huebner, 06/07/20
################################################################################

from glob import glob
import os
import re

import pysam
import pandas as pd

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/"

#### Libraries #################################################################
LIBRARIES = {'SP4210': ['D5864', 'F7638', 'F7639', 'F7640', 'F7641', 'F8851', 'F8852', 'F8853', 'F8854', 'F8855', 'F8856', 'F8857', 'F8858'],  # SP4210 - RBC1, Ritidian Beach Cave, Guam
             'SP4211': ['D5865', 'F7642', 'F7643', 'F7644', 'F7645', 'F8862', 'F8863', 'F8864', 'F8865', 'F8866', 'F8867'],  # SP4211 - RBC2, Ritidian Beach Cave, Guam
            }
LIBLIST = {lib: k for k, v in LIBRARIES.items() for lib in v}
CAPTURE_ARRAY = pd.read_csv("documentation/capture_arrays.tsv", sep="\t")
CAPTURE_RUNS = CAPTURE_ARRAY.loc[CAPTURE_ARRAY['array'].isin(['390K', '390Ksupp']), 'seqrun_id'].tolist()
################################################################################

#### Auxilliary functions ######################################################
def find_bams(sampleid, libid, seqrunids=CAPTURE_RUNS):
    """Find all BAM files per library and return list."""
    extract_seqrun = re.compile(r'SP421[01]-[A-Z][0-9]+-([0-9]+_[A-Z]+[0-9]+_[0-9]+_lane[0-9]).uniq.L35MQ25.bam')
    runids = [os.path.basename(fn).replace(".uniq.L35MQ25.bam", "")
              for fn in glob(f"analysis/{sampleid}/{sampleid}-{libid}*uniq.L35MQ25.bam")
              if extract_seqrun.search(os.path.basename(fn)).group(1) in seqrunids]  
    return runids
################################################################################

wildcard_constraints:
    lib = "[A-Z][0-9]+",
    type = "(all|deam)"

rule all:
    input: 
        "publication/supp_tables/SUM_nuclearcapture.csv"

#### Damageprofiler ############################################################

rule generate_hg19_bed:
    output:
        temp("analysis/tmp/hg19.bed")
    message: "Create BED file filtering for the human chromosomes"
    params:
        dict = "/mnt/solexa/Genomes/hg19_evan/whole_genome.dict"
    run:
        with open(output[0], "wt") as outfile:
            with open(params.dict, "rt") as dictfile:
                for i, line in enumerate(dictfile):
                    if line.startswith("@SQ") and i < 26:
                        _, chrom, length, _, _ = line.rstrip().split("\t")
                        outfile.write("{}\t0\t{}\n".format(chrom[3:], length[3:]))
                    outfile.write("gi|251831106|ref|NC_012920.1|\t0\t16569\n")

rule merge_bams_by_lib:
    input:
        "analysis/tmp/hg19.bed"
    output:
        "tmp/perlib_bams/{lib}.{type}.bam"
    message: "Merge all 1240K arrays of library {wildcards.lib}"
    params:
        sample = lambda wildcards: LIBLIST[wildcards.lib],
        bam = lambda wildcards: ",".join(find_bams(LIBLIST[wildcards.lib], wildcards.lib)),
        suffix = lambda wildcards: "uniq.L35MQ25.bam" if wildcards.type == "all" else "uniq.L35MQ25.deam.bam",
        reffasta = "/mnt/solexa/Genomes/hg19_evan/whole_genome.fa"
    shell:
        """
        samtools cat analysis/{params.sample}/{{{params.bam}}}.{params.suffix} | \
        samtools sort - | \
        samtools view -bh -L {input} - | \
        samtools calmd -u - {params.reffasta} > {output}
        """

rule remove_ambiguous_MD:
    # Fix problem of HTSJDK that it counts ambiguity bases, e.g. M, differently
    # than HTSLIB regarding the number of mismatches and therefore throws an
    # error for an illegal CIGAR string
    # see https://github.com/Integrative-Transcriptomics/DamageProfiler/issues/31
    input:
        "tmp/perlib_bams/{lib}.{type}.bam"
    output:
        temp("tmp/perlib_bams/{lib}.{type}.mdfix.bam")
    message: "Remove reads whose MD tag contains an IUPAC ambiguity code: {wildcards.lib} for {wildcards.type} sequencing data"
    run:
        bamfile = pysam.AlignmentFile(input[0], "rb")
        with pysam.AlignmentFile(output[0], "wb", template=bamfile) as outfile:
            for read in bamfile:
                correct_md = True
                if read.has_tag("MD"):
                    if "M" in read.get_tag("MD") or "R" in read.get_tag("MD"):
                        correct_md = False
                if correct_md:
                    outfile.write(read)

rule damageProfiler:
    input:
        "tmp/perlib_bams/{lib}.{type}.mdfix.bam"
    output:
        "tmp/damageprofiler/{lib}.{type}/identity_histogram.pdf"
    message: "Run DamageProfiler on library {wildcards.lib} for {wildcards.type} sequencing data"
    params:
        dir = "tmp/damageprofiler/{lib}.{type}",
        reffasta = "/mnt/solexa/Genomes/hg19_evan/whole_genome.fa",
        tmpdir = "tmp"
    threads: 4
    shell:
        """
        damageprofiler -Xms512M -Xmx12G \
             -Djava.awt.headless=true \
             -i {input} \
             -o {params.tmpdir} \
             -r {params.reffasta} \
             --title "{wildcards.lib}-{wildcards.type}"
        mv {params.tmpdir}/{wildcards.lib}.{wildcards.type}.mdfix/* {params.dir}
        rmdir {params.tmpdir}/{wildcards.lib}.{wildcards.type}.mdfix/
        """

rule summarise_damageprofiler:
    input:
        expand("tmp/damageprofiler/{lib}.{type}/identity_histogram.pdf", lib=LIBLIST.keys(), type=['all', 'deam'])
    output:
        "results/QUAL_damageprofile_perLib.RData"
    message: "Summarise the misincorporation statistics per library"
    params:
        dir = "tmp/damageprofiler"
    script:
        "scripts/SMRY_nuclearData-summarise_damageprofiler.R"

################################################################################

#### Conditional substitution ##################################################

rule summarise_condsubst:
    output:
        "results/QUAL_condSubst_perLib.RData"
    message: "Summarise the conditional substitution results per library"
    params:
        freq = "analysis/qual/conditionalSubsts_summary.txt",
        subst = "analysis/logs/5p3p_substitution_summary.txt"
    script:
        "scripts/SMRY_nuclearData-summarise_condsubst.R"

################################################################################

#### Summary ###################################################################

rule merge_library_summary:
    input:
        damage = "results/QUAL_damageprofile_perLib.RData",
        condsubst = "results/QUAL_condSubst_perLib.RData"
    output:
        "publication/supp_tables/SUM_nuclearcapture.csv"
    message: "Merge the summary of the library stats on nuclear data across samples"
    params:
        libdata = "results/SUM_libraries.csv",
        xhet = "analysis/qual/Xheterogeneity/X_contamination.txt"
    script:
        "scripts/SMRY_nuclearData-merge_library_summary.R"

################################################################################
