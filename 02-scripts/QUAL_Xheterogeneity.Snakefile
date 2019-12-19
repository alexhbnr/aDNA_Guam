################################################################################
# Project: Ancient human remains from Guam
# Part: Quality
# Step: Estimate contamination based on heterogeneity on X chromosome of males
#
# I will merge all sequencing runs of a single library and run ANGSD's
# contamination estimate based on the heterogeneity of the X chromosome for
# males.
#
# Alex Huebner, 17/11/19
################################################################################

from glob import glob
import os
import re

import numpy as np
import pandas as pd

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/"

#### SAMPLES ###################################################################
SAMPLES, = glob_wildcards("analysis/tmp/bam_per_lib/{sample}.bam")
################################################################################

rule all:
    input:
        expand("analysis/qual/Xheterogeneity/{male}_angsdXcont.log", male=SAMPLES),
        "analysis/qual/Xheterogeneity/X_contamination.txt"

rule count_X:
    output:  
        "/mnt/scratch/alexh/tmp/ancFlores/Xcont_{male}/angsdput.icnts.gz"
    message: "Generate binary count file for {wildcards.male}"
    params: 
        tmpdir = "/mnt/scratch/alexh/tmp/ancFlores/Xcont_{male}",
        bam = "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/analysis/tmp/bam_per_lib/{male}.bam"
    shell:
        """
        cd {params.tmpdir}
        angsd -i {params.bam} \
              -r X:5000000-154900000 \
              -doCounts 1 \
              -iCounts 1 \
              -minMapQ 25 \
              -minQ 30
        """

rule infer_X_cont:
    input:  
        "/mnt/scratch/alexh/tmp/ancFlores/Xcont_{male}/angsdput.icnts.gz"
    output: 
        "analysis/qual/Xheterogeneity/{male}_angsdXcont.log"
    message: "Estimate contamination on X chrosome for {wildcards.male}"
    params:
        bam = "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/analysis/tmp/bam_per_lib/{male}.bam"
    shell:
        """
        set +e
        /home/alexander_huebner/miniconda3/pkgs/angsd-0.923-h3ef6ad9_0/bin/contamination \
                    -a {input} \
                    -h /home/alexander_huebner/opt/angsd/RES/HapMapChrX.gz \
                    &> {output}
        rm -r $(dirname {input})
        if [ ! "$?" -eq 0 ]; then
            touch {output}
        fi
        """

rule summarise_X_cont:
    input:
        expand("analysis/qual/Xheterogeneity/{male}_angsdXcont.log", male=SAMPLES)
    output: 
        "analysis/qual/Xheterogeneity/X_contamination.txt"
    message: "Summarise the ANGSD contamination based from the X chromosome of males"
    run:
        extract_ML =  re.compile("Method1: new_llh Version: MoM:-*[0-9]\.[0-9]+ SE\(MoM\):[0-9\.\-ena]+ ML:-*([0-9]\.[0-9]+) SE\(ML\):([0-9]\.[0-9]+e[-\+][0-9]+)")
        decompose_id = re.compile("(SP[0-9]+)-([A-Z][0-9]+)_([a-z]+)_angsdXcont.log")

        estimates = []
        for fn in glob("analysis/qual/Xheterogeneity/*_angsdXcont.log"):
            sample, lib, flt = decompose_id.search(os.path.basename(fn)).groups()
            point = np.nan
            se = np.nan
            for line in open(fn, "rt"):
                if line.startswith("Method1: new_llh Version:"):
                    point, se = extract_ML.search(line.rstrip()).groups()
            estimates.append((sample, lib, flt, point, se))

        pd.DataFrame(estimates, columns=['sample', 'library', 'filter', 'point estimate', 'standard error']).to_csv(output[0], sep="\t", index=False, na_rep="NA")
