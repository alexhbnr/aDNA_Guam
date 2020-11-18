################################################################################
# Project: Ancient human remains from Guam
# Part: Data preparation
# Step: Upload sequencing data to ENA
#
# For uploading the sequencing data to ENA, I will go back to the unfiltered
# sequencing data after adapter-trimming by leeHOM and will demultiplex the
# files based on a perfect index match. I will merge the sequencing runs of
# same type (shotgun, MT capture, autosomal capture) for each library and fix the
# SAM header to only contain a single library.
#
# Alex Huebner, 08/11/2020
################################################################################

from glob import glob
import os

import pandas as pd

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam"

#### Sequencing runs ###########################################################
SEQRUNS = {'160818_SN7001204_0542_lane2': '/mnt/ngs_data/160818_SN7001204_0542_BHWCVKBCXX_R_PEdi_D5788_D1097_D5886/Bustard/Raw_Sequences/s_2_sequence.bam', 
           '160906_M02279_0022_lane1': '/mnt/ngs_data/160906_M02279_0022_000000000-AUCAV_BN_D2869/Bustard/Raw_Sequences/s_1_sequence.bam',
           '161208_SN7001204_0563_lane1': '/mnt/ngs_data/161208_SN7001204_0563_AH3WJVBCXY_R_PEdi_D4262_D4264/Bustard/Raw_Sequences/s_1_sequence.bam',
           '161208_SN7001204_0564_lane1': '/mnt/ngs_data/161208_SN7001204_0564_BH3WKLBCXY_R_PEdi_D4267_D4268/Bustard/Raw_Sequences/s_1_sequence.bam',
           '170630_D00829_0053_lane1': '/mnt/ngs_data/170630_D00829_0053_BHJ2FTBCXY_R_PEdi_F7703_F7509/Bustard/Raw_Sequences/s_1_sequence.bam',
           '170822_M02279_0141_lane1': '/mnt/ngs_data/170822_M02279_0141_000000000-BB4NL_BN_F5105/Bustard/Raw_Sequences/s_1_sequence.bam',
           '170906_D00829_0070_lane1': '/mnt/ngs_data/170906_D00829_0070_AHNVFNBCXY_R_PEdi_F8888_G1774/Bustard/Raw_Sequences/s_1_sequence.bam',
           '171020_D00829_0083_lane1': '/mnt/ngs_data/171020_D00829_0083_AHWTHYBCXY_R_PEdi_F5903_F5897/Bustard/Raw_Sequences/s_1_sequence.bam',
           '171020_D00829_0084_lane1': '/mnt/ngs_data/171020_D00829_0084_BHWVGYBCXY_R_PEdi_F5904_G3070/Bustard/Raw_Sequences/s_1_sequence.bam',
           '171110_D00829_0090_lane1': '/mnt/ngs_data/171110_D00829_0090_AH2GCVBCX2_R_PEdi_F5903_F5904/Bustard/Raw_Sequences/s_1_sequence.bam',
           '171110_D00829_0090_lane2': '/mnt/ngs_data/171110_D00829_0090_AH2GCVBCX2_R_PEdi_F5903_F5904/Bustard/Raw_Sequences/s_2_sequence.bam',
           '171110_D00829_0091_lane1': '/mnt/ngs_data/171110_D00829_0091_BH2M2TBCX2_R_PEdi_F5904/Bustard/Raw_Sequences/s_1_sequence.bam',
           '171110_D00829_0091_lane2': '/mnt/ngs_data/171110_D00829_0091_BH2M2TBCX2_R_PEdi_F5904/Bustard/Raw_Sequences/s_2_sequence.bam',
           }
ENRICHMENT = pd.read_csv("documentation/capture_arrays.tsv", sep="\t", index_col=[0])
ENRICHMENT['array'] = ENRICHMENT['array'].str.replace(r"390K(supp)*", "nuclear")
################################################################################

#### Libraries #################################################################
LIBRARIES = {'RBC1': ['D5864', 'F7638', 'F7639', 'F7640', 'F7641', 'F8851', 'F8852', 'F8853', 'F8854', 'F8855', 'F8856', 'F8857', 'F8858', 'F5861', 'F5056', 'F5057', 'F5058', 'F5059'],  # SP4210 - RBC1, Ritidian Beach Cave, Guam
             'RBC2': ['D5865', 'F7642', 'F7643', 'F7644', 'F7645', 'F8862', 'F8863', 'F8864', 'F8865', 'F8866', 'F8867', 'F5060', 'F5061', 'F5062', 'F5063'],  # SP4211 - RBC2, Ritidian Beach Cave, Guam
            }
# Generate list of all read groups in BAM files that are to be expected
LIBDICT = {lib: k for k, v in LIBRARIES.items() for lib in v}
LIBLIST = [lib for _, v in LIBRARIES.items() for lib in v]
################################################################################

#### Adapter removal without collapse and demultiplexing by perfect index match

localrules: index_list

rule demultiplexing:
    input: expand("data/splitBAM_raw/{seqrun_id}.done", seqrun_id=SEQRUNS.keys())

rule index_list:
    output:  
        "data/splitBAM_raw/{seqrun_id}.indexlist.txt"
    params: indexlist = lambda wildcards: os.path.dirname(SEQRUNS[wildcards.seqrun_id]) + "/../../build/lane" + wildcards.seqrun_id[-1] + "/indices.txt"
    run:
        indexlist = pd.read_csv(params.indexlist, sep="\t", index_col=[2])
        # Remove all index combinations that are not from Flores samples
        indexlist = indexlist.loc[indexlist.index.isin(LIBLIST)].dropna(axis=0, how="all")
        indexlist = indexlist.reset_index()
        indexlist.columns = ['readgroup', 'p7', 'p5']
        indexlist.to_csv(output[0], sep="\t", index=False)

rule leeHom_nocollapse:
    output:
        temp("data/splitBAM_raw/{seqrun_id}.nsorted.bam")
    message: "Run leeHOM on sequencing run {wildcards.seqrun_id} without collapsing partial overlapping samples"
    params:
        bam = lambda wildcards: SEQRUNS[wildcards.seqrun_id]
    log:
        leehom = "data/splitBAM_raw/{seqrun_id}.leeHOM.log",
        # filterReads = "data/splitBAM_raw/{seqrun_id}.filterReads.log"
    shell:
        """
        /mnt/solexa/bin/aLib//leehom/src/leeHom \
            -k ',' \
            -f 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTG' \
            -s 'GGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT' \
            -c 'ACACTCTTTCCCTACACGTCTGAACTCCAG,ACACTCTTTCCCACACGTCTGAACTCCAGT,ACACTCTTTCCCTACACACGTCTGAACTCC,CTCTTTCCCTACACGTCTGAACTCCAGTCA,GAAGAGCACACGTCTGAACTCCAGTCACII,GAGCACACGTCTGAACTCCAGTCACIIIII,GATCGGAAGAGCACACGTCTGAACTCCAGT,AGATCGGAAGAGCACACGTCTGAACTCCAG,AGAGCACACGTCTGAACTCCAGTCACIIII,ACACGTCTGAACTCCAGTCACIIIIIIIAT,GTGCACACGTCTGAACTCCAGTCACIIIII,AGCACACGTCTGAACTCCAGTCACIIIIII,CGTATGCCGTCTTCTGCTTGAAAAAAAAAA' \
            --log {log.leehom} \
            -u \
            -o /dev/stdout {params.bam} | \
        samtools sort -n -o {output} /dev/stdin
        """
        # /mnt/solexa/bin/aLib//pipeline/filterReads \
            # -o /dev/stdout \
            # --log {log.filterReads} /dev/stdin | \

rule split_BAM:
    input:
        indexlist="data/splitBAM_raw/{seqrun_id}.indexlist.txt",
        bam="data/splitBAM_raw/{seqrun_id}.nsorted.bam"
    output:
        touch("data/splitBAM_raw/{seqrun_id}.done")
    message: "Split BAM based on perfect index matching: {wildcards.seqrun_id}"
    params:
        dir = "data/splitBAM_raw/{seqrun_id}",
        splitBAM = "/mnt/genotyping/sk_pipelines/projects/aDNA_Flores/scripts/splitBAM.py"
    log:
        "data/splitBAM_raw/{seqrun_id}.splitBAM_raw.log"
    shell:
        """
        mkdir -p {params.dir}
        python {params.splitBAM} \
                -i {input.bam} \
                -l {input.indexlist} \
                --unaligned \
                -o {params.dir} > {log}
        """

################################################################################

#### Merge data belonging to the same library and enrichment type ##############

def evaluate_no_files(wildcards):
    """Evaluates the number of files that have to be generated."""
    filelist = pd.read_csv(checkpoints.create_folder_list.get(**wildcards).output[0], sep="\t")
    return [f"data/ena/{fn}.bam" for fn in filelist['library_type']]

def nofiles(libent):
    if os.path.isfile("data/splitBAM_raw/library_seqfiles.txt"):
        n = pd.read_csv("data/splitBAM_raw/library_seqfiles.txt", sep="\t", index_col=[0]) \
                .loc[libent]['filepaths'].count(",") + 1
    else:
        n = 0
    return n

def listfiles(libent):
    if os.path.isfile("data/splitBAM_raw/library_seqfiles.txt"):
        fns = pd.read_csv("data/splitBAM_raw/library_seqfiles.txt", sep="\t", index_col=[0]) \
                  .loc[libent]['filepaths'].replace(",", " ")
    else:
        fns = ""
    return fns

wildcard_constraints:
    libent = "[A-Z][0-9]+_[A-Za-z]+"

rule merge_enrichment_types:
    input:
        "data/ena/md5sum.txt"

checkpoint create_folder_list:
    input:
        expand("data/splitBAM_raw/{seqrun_id}.done", seqrun_id=SEQRUNS.keys())
    output:
        "data/splitBAM_raw/library_seqfiles.txt"
    message: "Create list of sequence files per library and enrichment type"
    params:
        dir = "data/splitBAM_raw"
    run:
        with open(output[0], "wt") as outfile:
            print("library_type\tfilepaths", file=outfile)
            files = {}
            for lib in LIBLIST:
                seqfiles = glob(f"{params.dir}/**/{lib}.bam", recursive=True)
                foldertypes = {'shotgun': [],
                            'MT': [],
                            'nuclear': []}
                for sfile in seqfiles:
                    foldertypes[ENRICHMENT.loc[os.path.basename(os.path.dirname(sfile))].iloc[0]].append(sfile)
                for enrichment in ['shotgun', 'MT', 'nuclear']:
                    if len(foldertypes[enrichment]) > 0:
                        print(f"{lib}_{enrichment}\t{','.join(foldertypes[enrichment])}", file=outfile)

rule concatenate_bams:
    input:
        "data/splitBAM_raw/library_seqfiles.txt"
    output:
        "data/ena/{libent}.raw.bam"
    message: "Concatenate BAM files of {wildcards.libent}"
    params:
        dir = "data/ena",
        nfiles = lambda wildcards: nofiles(wildcards.libent),
        fns = lambda wildcards: listfiles(wildcards.libent)
    shell:
        """
        mkdir -p {params.dir}
        if [[ {params.nfiles} -gt 1 ]]; then
            samtools cat -o {output} {params.fns}
        else
            ln -s ${{PWD}}/{params.fns} {output}
        fi
        """

rule fix_filter_BAM:
    input:
        "data/ena/{libent}.raw.bam"
    output:
        "data/ena/{libent}.flt.bam"
    message: "Fix the BAM file of sample {wildcards.libent} and remove quality failed reads"
    shell:
        """
        bam-fixpair --quiet {input} | \
        samtools view -bh -F 512 - > {output}
        """

rule add_rg_tag:
    input:
        "data/ena/{libent}.flt.bam"
    output:
        "data/ena/{libent}.bam"
    message: "Add RG tag for {wildcards.libent}"
    params:
        rg = lambda wildcards: f'-r ID:{wildcards.libent.split("_")[0]} -r LB:{wildcards.libent} -r SM:{LIBDICT[wildcards.libent.split("_")[0]]}'
    shell:
        """
        samtools addreplacerg {params.rg} {input} > {output}
        """

rule md5sum:
    input:
        evaluate_no_files
    output:
        "data/ena/md5sum.txt"
    message: "Calculate MD5sum of final BAM files"
    shell:
        "md5sum {input} > {output}"


################################################################################

#### Prepare CSV file for bulk upload ##########################################

rule sample_tables:
    input: 
        f"{workflow.basedir}/../05-results/ENA_submission_samples.tsv",
        # f"{workflow.basedir}/../05-results/ENA_submission_bam.tsv"

rule prepare_sample_table:
    output:
        f"{workflow.basedir}/../05-results/ENA_submission_samples.tsv"
    message: "Prepare bulk submission file for sampling registering at ENA"
    params: 
        template = f"{workflow.basedir}/../01-documentation/ENA_samples_template.tsv",
    run:
        # Read first three lines of template
        template = [line.rstrip() for line in open(params.template, "rt")]

        # Generate sample information list from Supplementary tables 1, 2, 4
        samplelist = pd.DataFrame.from_dict({'sample_alias': ['RBC1', 'RBC2']})
        samplelist = samplelist.assign(tax_id = 9606)
        samplelist = samplelist.assign(scientific_name = "Homo sapiens")
        samplelist = samplelist.assign(common_name = "human")

        # Write to file
        with open(output[0], "wt") as outfile:
            outfile.write(template[0] + "\n")
            outfile.write(template[1] + "\n")
        samplelist[['sample_alias', 'tax_id', 'scientific_name', 'common_name']]. \
                to_csv(output[0], sep="\t", index=False, mode="a")
        with open(output[0], "at") as outfile:
            outfile.write(template[3] + "\n")
            outfile.write(template[4] + "\n")


rule prepare_bam_table:
    input:
        "data/ena/md5sum.txt"
    output:
        f"{workflow.basedir}/../05-results/ENA_submission_bam.tsv"
    message: "Prepare bulk submission file for BAM files at ENA"
    run:
        # Read MD5 files
        ERS = {'RBC1': 'ERS5341448',
               'RBC2': 'ERS5341449'}
        md5 = pd.read_csv(input[0],
                          sep="\s+", header=None, names=['file_md5', 'fn'])
        md5['library_name'] = md5['fn'].str.extract(r'.+/([A-Z][0-9]+_[A-Za-z]+).bam')
        md5 = md5.drop(['fn'], axis=1)
        md5['file_name'] = "Pugach_PNAS_Guam/" + md5['library_name'] + ".bam"
        md5['libid'] = md5['library_name'].str.extract(r'([A-Z][0-9]+)_[A-Za-z]+')
        md5['type'] = md5['library_name'].str.extract(r'[A-Z][0-9]+_([A-Za-z]+)')
        md5['sample_alias'] = md5['libid'].replace(LIBDICT).replace(ERS)
        md5['library_construction_protocol'] = md5['type'].replace({'shotgun': "",
                                                                    'nuclear': "in-solution capture enrichment for human nuclear DNA",
                                                                    'MT': "in-solution capture enrichment for human mtDNA"})
        md5['library_strategy'] = md5['type'].replace({'shotgun': "WGS",
                                                       'nuclear': "Targeted-Capture",
                                                       'MT': "Targeted-Capture"})
        md5 = md5.assign(platform = "ILLUMINA")
        md5 = md5.assign(instrument_model = "Illumina HiSeq 2500")
        md5 = md5.assign(library_source = "GENOMIC")
        md5 = md5.assign(library_selection = "RANDOM")
        md5 = md5.assign(design_description = "")
        md5 = md5.assign(library_layout = "PAIRED")
        md5 = md5.assign(insert_size = 250)
        md5 = md5.assign()
        # Write to file
        md5[['sample_alias', 'platform', 'instrument_model', 'library_name', 'library_source',
             'library_selection', 'library_strategy',
             'design_description', 'library_construction_protocol',
             'library_layout', 'insert_size',
             'file_name', 'file_md5']]. \
                to_csv(output[0], sep="\t", index=False)


################################################################################
