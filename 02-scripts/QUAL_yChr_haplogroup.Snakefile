################################################################################
# Project: Ancient human remains from Guam
# Part: Quality
# Step: Y-chrosome haplogroup
#
# Alex Huebner, 17/11/19
################################################################################


workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/"

#### Samples ###################################################################
SAMPLES = ['SP4210']
################################################################################


rule all:
    input:
        "analysis/qual/yHaplo/haplogroups.ancGuam.Y.txt",
        "results/YCHR_hgvs.csv"

rule subset_vcf_to_Ychr:
    output:
        vcf = "analysis/qual/ancGuam.Y.vcf.gz",
        tbi = "analysis/qual/ancGuam.Y.vcf.gz.tbi",
    message: "Subset the random-allele sampled genotypes of sample SP4210 to the Y-chromosome sites"
    params:
        vcf = "analysis/genotypes/SP4210.all.vcf.gz"
    shell:
        """
        bcftools view -r Y {params.vcf} | \
        bcftools reheader -s SP4120.all | \
        bgzip > {output.vcf} 
        tabix {output.vcf}
        """

rule yhaplo:
    input:  
        "analysis/qual/ancGuam.Y.vcf.gz"
    output: 
        "analysis/qual/yHaplo/haplogroups.ancGuam.Y.txt"
    message: "Determine Y chromosome haplogroups using yhaplo"
    params:
        dir = "analysis/qual/yHaplo"
    shell:
        """
        set +u
        source /home/alexander_huebner/miniconda3/etc/profile.d/conda.sh
        conda activate py27
        set -u
        python /home/alexander_huebner/github/yhaplo/callHaplogroups.py \
                -i {input} \
                --ancStopThresh 1e6 \
                -asd -dsd -hpd \
                -o {params.dir}
        """

rule summarise_hgvs:
    input:
        "analysis/qual/yHaplo/haplogroups.ancGuam.Y.txt"
    output:
        "results/YCHR_hgvs.csv"
    message: "Summarise the output of yHaplo in HGVS format"
    params:
        dir = "analysis/qual/yHaplo"
    script:
        "scripts/QUAL_yChr_haplogroup-summarise_hgvs.R"
