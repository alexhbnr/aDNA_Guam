################################################################################
# Project: Ancient human remains from Guam
# Part: Summary
# Step: Summary of the major information on the number of SNPs per sample
#
# Alex Huebner, 18/11/19
################################################################################

from snakemake.utils import R
import pandas as pd

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/"

rule all:
    input:
        "results/SUM_SNPs.csv"

rule overlap_390K:
    output:
        "analysis/genotypes/dataset_info/390K.txt"
    message: "Summarise the number of SNPs per sample that overlap with the 390K array"
    params:
        array = "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/documentation/390K.bed.gz",
        vcf_dir = "analysis/genotypes"
    shell:
        """
        echo -e "sample\tfilter\tnSNPs" > {output}
        for vcf in {params.vcf_dir}/*.vcf.gz; do
            if [[ ${{vcf}} == *"deam"* ]]; then
                echo -e "$(basename ${{vcf}} | cut -d'.' -f1)\tdeam\t$(bcftools view -H -T {params.array} ${{vcf}} | wc -l)"
            else
                echo -e "$(basename ${{vcf}} | cut -d'.' -f1)\tall\t$(bcftools view -H -T {params.array} ${{vcf}} | wc -l)"
            fi
        done >> {output}
        """

rule overlap_390Ksupplement:
    output:
        "analysis/genotypes/dataset_info/390Ksupplement.txt"
    message: "Summarise the number of SNPs per sample that overlap with the 390K supplement array"
    params:
        array = "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/documentation/840K.bed.gz",
        vcf_dir = "analysis/genotypes"
    shell:
        """
        echo -e "sample\tfilter\tnSNPs" > {output}
        for vcf in {params.vcf_dir}/*.vcf.gz; do
            if [[ ${{vcf}} == *"deam"* ]]; then
                echo -e "$(basename ${{vcf}} | cut -d'.' -f1)\tdeam\t$(bcftools view -H -T {params.array} ${{vcf}} | wc -l)"
            else
                echo -e "$(basename ${{vcf}} | cut -d'.' -f1)\tall\t$(bcftools view -H -T {params.array} ${{vcf}} | wc -l)"
            fi
        done >> {output}
        """


rule overlap_archaicadmixture:
    output:
        "analysis/genotypes/dataset_info/archaicadmixture.txt"
    message: "Summarise the number of SNPs per sample that overlap with the archaic admxiture array"
    params:
        array = "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/documentation/Archaic.align.noN.sorted.bed.gz",
        vcf_dir = "analysis/genotypes"
    shell:
        """
        echo -e "sample\tfilter\tnSNPs" > {output}
        for vcf in {params.vcf_dir}/*.vcf.gz; do
            if [[ ${{vcf}} == *"deam"* ]]; then
                echo -e "$(basename ${{vcf}} | cut -d'.' -f1)\tdeam\t$(bcftools view -H -T {params.array} ${{vcf}} | wc -l)"
            else
                echo -e "$(basename ${{vcf}} | cut -d'.' -f1)\tall\t$(bcftools view -H -T {params.array} ${{vcf}} | wc -l)"
            fi
        done >> {output}
        """

rule overlap_HumanOrigins:
    output:
        "analysis/genotypes/dataset_info/HumanOrigins.txt"
    message: "Summarise the number of SNPs per sample that overlap with the Human Origins array"
    params:
        array = "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/documentation/humanorigins_snp.hg19.bed.gz",
        vcf_dir = "analysis/genotypes"
    shell:
        """
        echo -e "sample\tfilter\tnSNPs" > {output}
        for vcf in {params.vcf_dir}/*.vcf.gz; do
            if [[ ${{vcf}} == *"deam"* ]]; then
                echo -e "$(basename ${{vcf}} | cut -d'.' -f1)\tdeam\t$(bcftools view -H -T {params.array} ${{vcf}} | wc -l)"
            else
                echo -e "$(basename ${{vcf}} | cut -d'.' -f1)\tall\t$(bcftools view -H -T {params.array} ${{vcf}} | wc -l)"
            fi
        done >> {output}
        """

rule overlap_Affymetrix6:
    output:
        "analysis/genotypes/dataset_info/Affymetrix6.txt"
    message: "Summarise the number of SNPs per sample that overlap with the Affymetrix 6.0 array"
    params:
        array = "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/documentation/affimetrix6_snp.hg19.bed.gz",
        vcf_dir = "analysis/genotypes"
    shell:
        """
        echo -e "sample\tfilter\tnSNPs" > {output}
        for vcf in {params.vcf_dir}/*.vcf.gz; do
            if [[ ${{vcf}} == *"deam"* ]]; then
                echo -e "$(basename ${{vcf}} | cut -d'.' -f1)\tdeam\t$(bcftools view -H -T {params.array} ${{vcf}} | wc -l)"
            else
                echo -e "$(basename ${{vcf}} | cut -d'.' -f1)\tall\t$(bcftools view -H -T {params.array} ${{vcf}} | wc -l)"
            fi
        done >> {output}
        """

rule generate_BED_axiom_gw_ceu1:
    output:
        "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/documentation/Axiom_GW_CEU1.bed.gz"
    message: "Generate the BED file for the Affymetrix Axiom GWAS CEU1 array used by Hudjashov et al. 2017"
    params:
        annot = "/home/alexander_huebner/Documents/StonekingLab/390K_AffimetrixChips/Axiom_GW_CEU1.annotation.csv.gz"
    run:
        annot = pd.read_csv("/home/alexander_huebner/Documents/StonekingLab/390K_AffimetrixChips/Axiom_GW_CEU1.annotation.csv.gz",
                            sep=",",
                            comment="#",
                            usecols=[4, 5, 6])
        annot = annot.loc[annot['Physical Position'] != "---"]
        annot['Physical Position'] = annot['Physical Position'].astype(int)
        annot['Position End'] = annot['Position End'].astype(int)
        annot['Physical Position'] = annot['Physical Position'] - 1
        annot = annot.sort_values(['Chromosome', 'Physical Position'])
        annot.to_csv(output[0], sep="\t", index=False, header=False, compression="gzip")


rule overlap_axiom_gw_ceu1:
    input:
        "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/documentation/Axiom_GW_CEU1.bed.gz"
    output:
        "analysis/genotypes/dataset_info/Axiom_GW_CEU1.txt"
    message: "Summarise the number of SNPs per sample that overlap with the Affymetrix Axiom GW CEU1 array"
    params:
        vcf_dir = "analysis/genotypes"
    shell:
        """
        echo -e "sample\tfilter\tnSNPs" > {output}
        for vcf in {params.vcf_dir}/*.vcf.gz; do
            if [[ ${{vcf}} == *"deam"* ]]; then
                echo -e "$(basename ${{vcf}} | cut -d'.' -f1)\tdeam\t$(bcftools view -H -T {input} ${{vcf}} | wc -l)"
            else
                echo -e "$(basename ${{vcf}} | cut -d'.' -f1)\tall\t$(bcftools view -H -T {input} ${{vcf}} | wc -l)"
            fi
        done >> {output}
        """

rule generate_mega_bed:
    output:
        "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/documentation/MEGA.bed.gz"
    message: "Convert PLINK BIM file to BED for MEGA positions"
    params:
        bim = "/home/irina_pugach/aDNA_Indonesia_Guam/New_Guam/Data/MEGA/OGVP.Guam.bim"
    shell:
        """
        bioawk -t '{{print $1, $4 - 1, $4}}' {params.bim} | bgzip > {output}
        """

rule overlap_mega:
    input:
        "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/documentation/MEGA.bed.gz"
    output:
        "analysis/genotypes/dataset_info/MEGA.txt"
    message: "Summarise the number of SNPs per sample that overlap with the MEGA array"
    params:
        vcf_dir = "analysis/genotypes"
    shell:
        """
        echo -e "sample\tfilter\tnSNPs" > {output}
        for vcf in {params.vcf_dir}/*.vcf.gz; do
            if [[ ${{vcf}} == *"deam"* ]]; then
                echo -e "$(basename ${{vcf}} | cut -d'.' -f1)\tdeam\t$(bcftools view -H -T {input} ${{vcf}} | wc -l)"
            else
                echo -e "$(basename ${{vcf}} | cut -d'.' -f1)\tall\t$(bcftools view -H -T {input} ${{vcf}} | wc -l)"
            fi
        done >> {output}
        """

rule summary:
    input:
        ca_390K = "analysis/genotypes/dataset_info/390K.txt",
        ca_390Ksupp = "analysis/genotypes/dataset_info/390Ksupplement.txt",
        archaicadmix = "analysis/genotypes/dataset_info/archaicadmixture.txt",
        human_origins = "analysis/genotypes/dataset_info/HumanOrigins.txt",
        affymetrix6 = "analysis/genotypes/dataset_info/Affymetrix6.txt",
        axiom_gw_ceu1 = "analysis/genotypes/dataset_info/Axiom_GW_CEU1.txt",
        mega = "analysis/genotypes/dataset_info/MEGA.txt"
    output:
        "results/SUM_SNPs.csv"
    message: "Summarises all information regarding the libraries of the Flores samples"
    params: 
        dir = "analysis/genotypes/dataset_info"
    script:
        "scripts/SMRY_samples-summary.R"
