################################################################################
# Project: Ancient human remains from Guam
# Part: Data preparation
# Step: Filter sequencing data
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
BAMS = {os.path.basename(fn)[:-4]: sample
        for sample in SAMPLES
        for fn in glob(f"analysis/{sample}/*.bam") 
        if "L35MQ25" not in fn}
################################################################################

#### Parameters ################################################################
MINLEN = 35
MINMQ = 25
################################################################################

wildcard_constraints:
    sample = "[A-Z]+[0-9]+",
    id = "[A-Z]+[0-9]+-[A-Z][0-9]+-[A-Za-z0-9_]+"

rule all:
    input: ["analysis/{sample}/{id}.uniq.L35MQ25.deam.bam".format(sample=BAMS[bam], id=bam) for bam in BAMS.keys()],
           "analysis/logs/average_length.txt",
           "analysis/logs/5p3p_substitution_summary.txt",
           "analysis/logs/deam/noReads_deaminatedonly.txt"

rule analyseBAM:
    # Remove unmapped, non-merged, filter-flagged sequences and create summary statistic
    output: 
        temp("analysis/{sample}/{id}.L35MQ25.bam")
    message: "AnalyzeBAM: BAM file {wildcards.id} of sample {wildcards.sample}"
    params:
        analyzeBAM = "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/scripts/analyzeBAM.py",
        bam = "analysis/{sample}/{id}.bam",
        sumstats = lambda wildcards: "analysis/logs/sumstats/{}-{}-L35MQ25.txt".format(wildcards.sample, wildcards.id[7:])
    shell:
        """
        mkdir -p analysis/logs/sumstats
        python {params.analyzeBAM} \
                --input {params.bam} \
                --minlength {MINLEN} \
                --minqual {MINMQ} \
                --statreport {params.sumstats} \
                --sample {wildcards.id} \
                --output analysis/{wildcards.sample}/{wildcards.id}
        """

rule sort_length_mq_filtered_bam:
    input:
        "analysis/{sample}/{id}.L35MQ25.bam"
    output:
        temp("analysis/{sample}/{id}.L35MQ25.sorted.bam")
    message: "Sort BAM file filtered for minimal length and mapping quality: {wildcards.id}"
    shell:
        """
        samtools sort -o {output} {input}
        """

rule rmdup_dedup:
    input:
        "analysis/{sample}/{id}.L35MQ25.sorted.bam"
    output:
        "analysis/{sample}/{id}.uniq.L35MQ25.bam"
    message: "Remove duplicate reads using DeDup: {wildcards.id}"
    params:
        tmpdir = "analysis/tmp"
    shell:
        """
        mkdir -p analysis/logs/dedup
        mkdir -p {params.tmpdir}
        dedup -i {input} \
              --merged \
              -o {params.tmpdir}/
        mv {params.tmpdir}/{wildcards.id}.L35MQ25.sorted.hist analysis/logs/dedup/{wildcards.id}.L35MQ25.hist
        mv {params.tmpdir}/{wildcards.id}.L35MQ25.sorted.log analysis/logs/dedup/{wildcards.id}.L35MQ25.log
        samtools sort -o {output} {params.tmpdir}/{wildcards.id}.L35MQ25.sorted_rmdup.bam
        rm {params.tmpdir}/{wildcards.id}.L35MQ25.sorted_rmdup.bam {params.tmpdir}/{wildcards.id}.L35MQ25.sorted.dedup.json
        """

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

rule samtools_calmd:
    # Remove all decoy chromosomes and fix the MD field
    input:
        bam = "analysis/{sample}/{id}.uniq.L35MQ25.bam",
        bed = "analysis/tmp/hg19.bed"
    output:
        temp("analysis/tmp/{sample}_{id}.uniq.L35MQ25.bam")
    message: "Fix MD tags using samtools: {wildcards.id}"
    params:
        reffasta = "/mnt/solexa/Genomes/hg19_evan/whole_genome.fa"
    shell:
        """
        samtools view -bh -L {input.bed} {input.bam} | \
        samtools calmd -u - {params.reffasta} > {output}
        """

rule remove_ambiguous_MD:
    # Fix problem of HTSJDK that it counts ambiguity bases, e.g. M, differently
    # than HTSLIB regarding the number of mismatches and therefore throws an
    # error for an illegal CIGAR string
    # see https://github.com/Integrative-Transcriptomics/DamageProfiler/issues/31
    input:
        "analysis/tmp/{sample}_{id}.uniq.L35MQ25.bam"
    output:
        temp("analysis/tmp/{sample}_{id}.uniq.L35MQ25.mdfix.bam")
    message: "Remove reads whose MD tag contains an IUPAC ambiguity code: {wildcards.sample}"
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
        "analysis/tmp/{sample}_{id}.uniq.L35MQ25.mdfix.bam"
    output:
        "analysis/logs/damageprofiler/{sample}/{id}/identity_histogram.pdf"
    message: "Run DamageProfiler on library {wildcards.id} of sample {wildcards.sample}"
    params:
        dir = "analysis/logs/damageprofiler/{sample}/{id}",
        reffasta = "/mnt/solexa/Genomes/hg19_evan/whole_genome.fa",
        tmpdir = "analysis/tmp"
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
        mv {params.tmpdir}/{wildcards.sample}_{wildcards.id}.uniq.L35MQ25.mdfix/* {params.dir}
        rmdir {params.tmpdir}/{wildcards.sample}_{wildcards.id}.uniq.L35MQ25.mdfix/
        """

rule average_length:
    input: 
        ["analysis/logs/damageprofiler/{sample}/{id}/identity_histogram.pdf".format(sample=BAMS[bam], id=bam) for bam in BAMS.keys()]
    output:
        "analysis/logs/average_length.txt"
    message: "Summarise the read length distribution of all libraries"
    params:
        dir =  "analysis/logs/damageprofiler"
    run:
        R("""
          library(data.table)
          library(tidyverse)
          
          length_profiles <- map_df(list.files(path = "{params.dir}"),
                                    function(s) {{
                               map_df(list.files(path = paste("{params.dir}", 
                                                              s, sep = "/"), full.names = T), function(l) {{
                                 lgdist <- fread(paste0(l, "/lgdistribution.txt")) %>%
                                           group_by(Length) %>%
                                           summarise(Occurrences = sum(Occurrences))
                                 lengths <- rep(lgdist$Length, lgdist$Occurrences)
                                 tibble(`sample` = s,
                                        `library` = basename(l),
                                        mean = mean(lengths),
                                        median = median(lengths),
                                        mode = lgdist %>%
                                               arrange(desc(Occurrences)) %>%
                                               pull(Length) %>% .[1],
                                        max = max(lengths))
                               }})
                             }})
          
          fwrite(length_profiles, sep = "\t", file = "{output}")
        """)

rule summarise_ct_frequency:
    input: 
        ["analysis/logs/damageprofiler/{sample}/{id}/identity_histogram.pdf".format(sample=BAMS[bam], id=bam) for bam in BAMS.keys()]
    output:
        "analysis/logs/5p3p_substitution_summary.txt"
    message: "Summarise the C to T frequency at first ten bases of each the 5' and 3' end"
    params:
        dir = "analysis/logs/damageprofiler"
    run:
        R("""
          library(data.table)
          library(tidyverse)
          
          subst_profiles <- map_df(list.files(path = "{params.dir}"),
                                    function(s) {{
                               map_df(list.files(path = paste("{params.dir}", 
                                                              s, sep = "/"), full.names = T), function(l) {{
                                 fiveprime <- fread(paste0(l, "/5pCtoT_freq.txt"),
                                                    skip = 3) %>%
                                              slice(1:10) %>%
                                              mutate(`sample` = s,
                                                     `library` = basename(l)) %>%
                                              spread(pos, `5pC>T`)
                                 threeprime <- fread(paste0(l, "/3p_freq_misincorporations.txt"),
                                                     skip = 3, select = 1:2) %>%
                                               slice(16:25) %>%
                                               mutate(`sample` = s,
                                                      `library` = basename(l),
                                                      Pos = (-1) * (Pos + 1)) %>%
                                              spread(Pos, `C>T`)
                                 left_join(fiveprime, threeprime, by=c("sample" = "sample",
                                                                       "library"  = "library"))
                               }})
                             }})
          
          fwrite(subst_profiles, sep = "\t", file = "{output}")
        """)

rule filterBAM:
    # Filtering the non-UDG treated samples for C-T subsitutions on the
    # last three positions of a read
    input: 
        "analysis/{sample}/{id}.uniq.L35MQ25.bam"
    output: 
        "analysis/{sample}/{id}.uniq.L35MQ25.deam.bam"
    message: "Filter library {wildcards.id} of sample {wildcards.sample} for C-T substitutions"
    log: "analysis/logs/deam/{sample}_{id}-filterBAM.log"
    params:
        filterBAM = "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/scripts/filterBAM.py",
    shell:
        """
        python {params.filterBAM} \
                -p3 0,-1,-2 \
                -p5 0,1,2 \
                -suffix "deam" \
                {input} 2> {log}
        """

rule summarise_noReads:
    input:
        ["analysis/{sample}/{id}.uniq.L35MQ25.deam.bam".format(sample=BAMS[bam], id=bam) for bam in BAMS.keys()]
    output:
        "analysis/logs/deam/noReads_deaminatedonly.txt"
    message: "Summarise the number of reads with terminal C-to-T subsitutions per sample"
    shell:
        """
        echo -e "id\tnSamples" > {output}
        for bam in {input}; do
            id=$(echo $(basename $(echo ${{bam}} | cut -f 6-) .uniq.L35MQ25.deam.bam))
            echo -e "${{id}}\t$(samtools view -c ${{bam}})"
        done >> {output}
        """

################################################################################
