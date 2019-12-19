################################################################################
# Project: Ancient human remains from Guam
# Part: Data preparation
# Step: Annotate and filter genotypes for each capture array
#
# Genotypes were generated for all sites that were covered with at least one
# base. In the following, the genotypes of the samples are merged by read filter
# status (all, deam) and genotyping method (MPISHH, MPIEVA) and subset to the
# array of interest (1240K, archaic admixture). The observed genotypes at the
# sites targeted by the array are compared to the expected ones provided by the
# annotation file. Genotypes different from the expected ones are set to
# missing and the alternative alleles are fixed. Finally, the VCF file is
# converted into EIGENSTRAT format.
#
# Alex Huebner, 18/11/19
################################################################################

import allel

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam"

#### Samples ###################################################################
SAMPLES, = glob_wildcards("analysis/bam/{sample}.all.nuclear.bam")
################################################################################

#### Genotype info #############################################################
FLT = ['all', 'deam']
ARRAY = {'1240K': "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/documentation/390Kplus840K.bed.gz",
         'archaicAdmix': "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/documentation/Archaic.align.noN.sorted.bed.gz"}
SNPINFO = {'1240K': "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/documentation/1240K.snp.gz",
           'archaicAdmix': "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/documentation/archaicAdmix.snp.gz"}
################################################################################

wildcard_constraints:
    array = "(1240K|archaicAdmix)",
    sample = "SP[0-9]+(HQ)*",
    flt = "(all|deam)"

rule all:
    input: 
        expand("analysis/genotypes/{array}/ancGuam.{flt}.geno", array=ARRAY, flt=FLT)

rule subset_VCF_to_array_sites:
    output:
        vcf = temp("analysis/genotypes/{array}/{sample}.{flt}.vcf.gz"),
        tbi = temp("analysis/genotypes/{array}/{sample}.{flt}.vcf.gz.tbi")
    message: "Subset genotypes of sample {wildcards.sample} filtered for {wildcards.flt} reads to sites of {wildcards.array} array" 
    params: 
        vcf = "analysis/genotypes/{sample}.{flt}.vcf.gz",
        bed = lambda wildcards: ARRAY[wildcards.array]
    shell:
        """
        bcftools view -T {params.bed} {params.vcf} | \
        bcftools norm -d all -O z -o {output.vcf} -
        tabix {output.vcf}
        """

rule merge_VCF:
    input:
        vcf = lambda wildcards: [f"analysis/genotypes/{wildcards.array}/{sample}.{wildcards.flt}.vcf.gz" for sample in SAMPLES],
        tbi = lambda wildcards: [f"analysis/genotypes/{wildcards.array}/{sample}.{wildcards.flt}.vcf.gz.tbi" for sample in SAMPLES]
    output:
        "analysis/genotypes/{array}/ancGuam.{flt}.raw.vcf.gz"
    message: "Merge all samples for array {wildcards.array} using {wildcards.flt} reads"
    shell:
        """
        bcftools merge -m all -O z -o {output} {input.vcf}
        """

rule annotate:
    input:
        "analysis/genotypes/{array}/ancGuam.{flt}.raw.vcf.gz"
    output:
        temp("analysis/genotypes/{array}/ancGuam.{flt}.annot_nonsorted.vcf.gz")
    message: "Annotate VCF with meta-information of the {wildcards.array} array generated from {wildcards.flt} reads"
    params:
        annotateVCFwSNPfile = "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/scripts/annotateVCFwSNPfile.py",
        snp = lambda wildcards: SNPINFO[wildcards.array]
    shell:
        """
        {params.annotateVCFwSNPfile} \
            -i {input} \
            -o {output} \
            -a {params.snp}
        """

rule sort_vcf_file:
    input:
        "analysis/genotypes/{array}/ancGuam.{flt}.annot_nonsorted.vcf.gz"
    output:
        "analysis/genotypes/{array}/ancGuam.{flt}.annot.vcf.gz"
    message: "Sort VCF file of the {wildcards.array} array generated from {wildcards.flt} reads"
    shell:
        """
        bcftools sort -O z -o {output} {input}
        """

rule generate_map:
    input:
        "analysis/genotypes/{array}/ancGuam.{flt}.annot.vcf.gz"
    output:
        "analysis/genotypes/{array}/ind.{flt}.map"
    message: "Generate individual map for Eigenstrat conversion"
    run:
        vcf = allel.read_vcf(input[0], fields=['samples'])
        with open(output[0], "wt") as outfile:
            for sm in vcf['samples']:
                outfile.write(f"{sm}\tU\t{sm}\n")

rule convert_to_eigenstrat:
    input:
        vcf = "analysis/genotypes/{array}/ancGuam.{flt}.annot.vcf.gz",
        map = "analysis/genotypes/{array}/ind.{flt}.map"
    output:
        "analysis/genotypes/{array}/ancGuam.{flt}.geno"
    message: "Convert VCF captured on {wildcards.array} array generated from {wildcards.flt} reads"
    params:
        vcf2eigenstrat = "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/scripts/vcf2eigenstrat.jl",
        outputprefix = "analysis/genotypes/{array}/ancGuam.{flt}"
    shell:
        """
        julia {params.vcf2eigenstrat} \
            --input {input.vcf} \
            --output {params.outputprefix} \
            --map {input.map}
        """
