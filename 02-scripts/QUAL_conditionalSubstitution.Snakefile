################################################################################
# Project: Ancient human remains from Guam
# Part: Quality
# Step: Analysis for conditional substitution patterns
#
# The analysis evaluates the frequency of C-to-T substitutions at reads that
# have a C-to-T substitution at the terminal base at both read ends.
#
# Alex Huebner, 17/11/19
################################################################################

from snakemake.utils import R
from glob import glob
import os.path

import pysam

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam"

#### Samples ###################################################################
SAMPLES = ['SP4210', 'SP4211']
BAMS = {os.path.basename(fn).replace(".uniq.L35MQ25.bam", ""): sample
        for sample in SAMPLES
        for fn in glob(f"analysis/{sample}/*.uniq.L35MQ25.bam")
        if os.path.basename(fn).replace(".uniq.L35MQ25.bam", "") != "SP4210-F8854-170906_D00829_0070_lane1"}  # no conditional reads
################################################################################

def determine_reffasta(bamfn):
    ''' Determine the reference genome based on the SQ entries in the BAM header
    '''
    if os.path.isfile(bamfn):
        references = pysam.AlignmentFile(bamfn).header['SQ']
    else:
        references = []
    if len(references) > 1:
        return "/mnt/solexa/Genomes/hg19_evan/whole_genome.fa"
    else:
        return "/mnt/solexa/Genomes/human_MT/whole_genome.fa"


wildcard_constraints:
    sample = "[A-Z]+[0-9]+",
    id = "[A-Z]+[0-9]+-[A-Z][0-9]+-[A-Za-z0-9_]+"

localrules: link_bam, generate_hg19_bed

rule all:
    input: 
        "analysis/qual/conditionalSubsts_summary.txt",
        "analysis/qual/conditionalSubsts_contaminationEstimate.txt"

rule link_bam:
    output:
        "/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.bam"
    message: "Link the BAM {wildcards.id} from sample {wildcards.sample} to tmp"
    params:
        bam = "analysis/{sample}/{id}.uniq.L35MQ25.bam"
    shell:
        "ln -s ${{PWD}}/{params.bam} {output}"

rule generate_hg19_bed:
    output:
        temp("/mnt/scratch/alexh/tmp/ancGuam/contSubst/hg19.bed")
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

rule filterBAM_5p:
    # Filtering the non-UDG treated samples for C-T subsitutions on the
    # last position of the 5' end of a read
    input: 
        "/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.bam"
    output: 
        "/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.deam_5p.bam"
    message: "Filter library {wildcards.id} of sample {wildcards.sample} for C-T substitutions on the 5' end"
    log: "/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.deam_5p-filterBAM.log"
    params:
        filterBAM = "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/scripts/filterBAM.py",
    shell:
        """
        python {params.filterBAM} \
                -p5 0 \
                -suffix "deam_5p" \
                {input} >> {log}
        """

rule samtools_calmd_5p:
    # Remove all decoy chromosomes and fix the MD field
    input:
        bam = "/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.deam_5p.bam",
        bed = "/mnt/scratch/alexh/tmp/ancGuam/contSubst/hg19.bed"
    output:
        temp("/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.deam_5p.calmd.bam")
    message: "Fix MD tags using samtools for the 5' end: {wildcards.id}"
    params:
        reffasta = lambda wildcards: determine_reffasta(f"/mnt/scratch/alexh/tmp/ancGuam/contSubst/{wildcards.sample}/{wildcards.id}.deam_5p.bam")
    shell:
        """
        samtools view -bh -L {input.bed} {input.bam} | \
        samtools calmd -u - {params.reffasta} > {output}
        """

rule remove_ambiguous_MD_5p:
    # Fix problem of HTSJDK that it counts ambiguity bases, e.g. M, differently
    # than HTSLIB regarding the number of mismatches and therefore throws an
    # error for an illegal CIGAR string
    # see https://github.com/Integrative-Transcriptomics/DamageProfiler/issues/31
    input:
        "/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.deam_5p.calmd.bam"
    output:
        temp("/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.deam_5p.calmd.mdfix.bam")
    message: "Remove reads whose MD tag contains an IUPAC ambiguity code: {wildcards.sample} of library {wildcards.id}"
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


rule damageProfiler_5p:
    input:
        "/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.deam_5p.calmd.mdfix.bam"
    output:
        "/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.deam_5p/identity_histogram.pdf"
    message: "Run DamageProfiler on library {wildcards.id} of sample {wildcards.sample} filter for C-to-T substitutions on the 5' end"
    params:
        dir = "/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}",
        reffasta = lambda wildcards: determine_reffasta(f"/mnt/scratch/alexh/tmp/ancGuam/contSubst/{wildcards.sample}/{wildcards.id}.deam_5p.bam"),
        tmpdir = "/mnt/scratch/alexh/tmp/ancGuam/contSubst/tmp"
    threads: 4
    shell:
        """
        java -Xms512M -Xmx12G \
             -Djava.awt.headless=true \
             -jar /home/alexander_huebner/miniconda3/share/damageprofiler-0.4.5-1/DamageProfiler-0.4.5.jar \
             -i {input} \
             -o {params.tmpdir} \
             -r {params.reffasta} \
             --title {wildcards.id}
        mv {params.tmpdir}/{wildcards.id}.deam_5p.calmd.mdfix/* {params.dir}/{wildcards.id}.deam_5p/
        rmdir {params.tmpdir}/{wildcards.id}.deam_5p.calmd.mdfix/
        """

rule filterBAM_3p:
    # Filtering the non-UDG treated samples for C-T subsitutions on the
    # last position of the 3' end of a read
    input: 
        "/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.bam"
    output: 
        "/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.deam_3p.bam"
    message: "Filter library {wildcards.id} of sample {wildcards.sample} for C-T substitutions on the 3' end"
    log: "/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.deam_3p-filterBAM.log"
    params:
        filterBAM = "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/scripts/filterBAM.py",
    shell:
        """
        python {params.filterBAM} \
                -p3 0 \
                -suffix "deam_3p" \
                {input} >> {log}
        """

rule samtools_calmd_3p:
    # Remove all decoy chromosomes and fix the MD field
    input:
        bam = "/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.deam_3p.bam",
        bed = "/mnt/scratch/alexh/tmp/ancGuam/contSubst/hg19.bed"
    output:
        temp("/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.deam_3p.calmd.bam")
    message: "Fix MD tags using samtools for the 3' end: {wildcards.id}"
    params:
        reffasta = lambda wildcards: determine_reffasta(f"/mnt/scratch/alexh/tmp/ancGuam/contSubst/{wildcards.sample}/{wildcards.id}.deam_3p.bam")
    shell:
        """
        samtools view -bh -L {input.bed} {input.bam} | \
        samtools calmd -u - {params.reffasta} > {output}
        """

rule remove_ambiguous_MD_3p:
    # Fix problem of HTSJDK that it counts ambiguity bases, e.g. M, differently
    # than HTSLIB regarding the number of mismatches and therefore throws an
    # error for an illegal CIGAR string
    # see https://github.com/Integrative-Transcriptomics/DamageProfiler/issues/31
    input:
        "/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.deam_3p.calmd.bam"
    output:
        temp("/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.deam_3p.calmd.mdfix.bam")
    message: "Remove reads whose MD tag contains an IUPAC ambiguity code: {wildcards.sample} of library {wildcards.id}"
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

rule damageProfiler_3p:
    input:
        "/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.deam_3p.calmd.mdfix.bam"
    output:
        "/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.deam_3p/identity_histogram.pdf"
    message: "Run DamageProfiler on library {wildcards.id} of sample {wildcards.sample} filter for C-to-T substitutions on the 3' end"
    params:
        dir = "/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}",
        reffasta = lambda wildcards: determine_reffasta(f"/mnt/scratch/alexh/tmp/ancGuam/contSubst/{wildcards.sample}/{wildcards.id}.deam_3p.bam"),
        tmpdir = "/mnt/scratch/alexh/tmp/ancGuam/contSubst/tmp"
    threads: 4
    shell:
        """
        java -Xms512M -Xmx12G \
             -Djava.awt.headless=true \
             -jar /home/alexander_huebner/miniconda3/share/damageprofiler-0.4.5-1/DamageProfiler-0.4.5.jar \
             -i {input} \
             -o {params.tmpdir} \
             -r {params.reffasta} \
             --title {wildcards.id}
        mv {params.tmpdir}/{wildcards.id}.deam_3p.calmd.mdfix/* {params.dir}/{wildcards.id}.deam_3p/
        rmdir {params.tmpdir}/{wildcards.id}.deam_3p.calmd.mdfix/
        """

rule summarise_condsubst:
    input:
        p3 = ["/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.deam_5p/identity_histogram.pdf".format(sample=BAMS[bam], id=bam) for bam in BAMS.keys()],
        p5 = ["/mnt/scratch/alexh/tmp/ancGuam/contSubst/{sample}/{id}.deam_3p/identity_histogram.pdf".format(sample=BAMS[bam], id=bam) for bam in BAMS.keys()]
    output:
        "analysis/qual/conditionalSubsts_summary.txt"
    message: "Summarise the conditional substitution patterns"
    params:
        dir = "/mnt/scratch/alexh/tmp/ancGuam/contSubst"
    run:
        R("""
          library(data.table)
          library(tidyverse)

          samples <- list.files(path = "{params.dir}", full.names = T)
          libraries <- sapply(samples, function (s) {{
                         lib_paths <- str_replace_all(basename(list.files(path = s,
                                                                          pattern = ".+\\\.deam_3p$")),
                                                      "\\\.deam_3p", "")
                       }}, USE.NAMES = F) %>%
                       unlist()
          
          conditionalsubst_profiles <- map_df(libraries, function(l) {{
              # Reads with C-T substitution on 3' end
              threeprime <- fread(paste0("{params.dir}", "/",
                                         str_sub(l, 1, 6), "/", l,
                                         ".deam_3p/5p_freq_misincorporations.txt"),
                                  skip = 3, select = 1:2) %>%
                                  mutate(`sample` = str_sub(l, 1, 6),
                                         `library` = l,
                                         Pos  = Pos + 1) %>%
                                  slice(1:2) %>%
                                  spread(Pos, `C>T`)
              # Reads with C-T substitution on 5' end
              fiveprime <- fread(paste0("{params.dir}", "/",
                                         str_sub(l, 1, 6), "/", l,
                                         ".deam_5p/3p_freq_misincorporations.txt"),
                                  skip = 3, select = 1:2) %>%
                                  mutate(`sample` = str_sub(l, 1, 6),
                                         `library` = l,
                                         Pos  = (-1) * (Pos + 1)) %>%
                                  slice(24:25) %>%
                                  spread(Pos, `C>T`)
              left_join(threeprime, fiveprime,
                        by = c("sample" = "sample",
                               "library" = "library"))
            }})
          
          fwrite(conditionalsubst_profiles, sep = "\t", file = "{output}")
        """)

rule calculate_contamination:
    input:
        "analysis/qual/conditionalSubsts_summary.txt"
    output:
        "analysis/qual/conditionalSubsts_contaminationEstimate.txt"
    message: "Estimate the contamination based on observed conditional substitution rate compared to the overall substitution rate"
    params:
        subst = "analysis/logs/5p3p_substitution_summary.txt"
    run:
        R("""
          library(data.table)
          library(tidyverse)

          condsubst <- fread("{input}",
                             header = T,
                             select = c(1, 2, 3, 6),
                             col.names = c("sample", "library", "cond5p", "cond3p"))
          allsubst <- fread("{params.subst}",
                             header = T,
                             select = c(1, 2, 3, 22),
                             col.names = c("sample", "library", "all5p", "all3p"))

          contamination <- left_join(allsubst, condsubst,
                                     by = c("sample" = "sample",
                                            "library" = "library")) %>%
                           mutate(cont3p = 1 - (all3p / cond3p),
                                  cont5p = 1 - (all5p / cond5p)) %>%
                           mutate(contamination = apply(.[, c(7, 8)], 1, mean)) %>%
                           select(sample, library, contamination)
          fwrite(contamination, sep = "\t", file = "{output}")
        """)
