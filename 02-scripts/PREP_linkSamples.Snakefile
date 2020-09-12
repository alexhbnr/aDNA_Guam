################################################################################
# Project: Ancient human remains from Guam
# Part: Data preparation
# Step: Link de-multiplexed BAM files into sequence folders
#
# Alex Huebner, 17/11/19
################################################################################

from glob import glob
import os.path
import re
import sys

workdir: "/mnt/genotyping/sk_pipelines/projects/aDNA_Guam"

#### Libraries #################################################################
LIBRARIES = {'SP4210': ['D5864', 'F7638', 'F7639', 'F7640', 'F7641', 'F8851', 'F8852', 'F8853', 'F8854', 'F8855', 'F8856', 'F8857', 'F8858', 'F5861', 'F5056', 'F5057', 'F5058', 'F5059'],  # SP4210 - RBC1, Ritidian Beach Cave, Guam
             'SP4211': ['D5865', 'F7642', 'F7643', 'F7644', 'F7645', 'F8862', 'F8863', 'F8864', 'F8865', 'F8866', 'F8867', 'F5060', 'F5061', 'F5062', 'F5063'],  # SP4211 - RBC2, Ritidian Beach Cave, Guam
            }
LIBLIST = {lib: k for k, v in LIBRARIES.items() for lib in v}
################################################################################

#### File list #################################################################
def determine_filename(lib, fns):
    ''' Return adapted filename of scheme <sample>-<library>-<seqrun>.bam
    '''
    sample = LIBLIST[lib]
    seqruns = [re.search("data/splitBAM/([0-9]+_[A-Z]+[0-9]+_[0-9]+_lane[0-9])/[A-Z][0-9]+.bam", fn).group(1)
               for fn in fns]
    return [f"analysis/{sample}/{sample}-{lib}-{seqrun}.bam" for seqrun in seqruns]
# Identify BAM generated by de-multiplexing
FILELIST = [glob(f"data/splitBAM/**/{lib}.bam") for lib in LIBLIST]
# Generate list of new filenames based on found BAMs
EXPECTED_FILELIST = [determine_filename(lib, glob(f"data/splitBAM/**/{lib}.bam"))
                     for lib in LIBLIST]
# Generate dict with new and old filenames
FILEDICT = {efnlist[i]: fnlist[i]
            for fnlist, efnlist in zip(FILELIST, EXPECTED_FILELIST)
            for i in range(0, len(fnlist))}
    

rule link_bamfiles:
    output:
        [bamfn for bamfn in FILEDICT.keys()]
    message: "Link BAM files"
    run:
        for fn in FILEDICT:
            if not os.path.isfile(fn):
                print("ln -s ${{PWD}}/{} {}".format(FILEDICT[fn], fn),
                    file=sys.stderr)
                subprocess.run("ln -s ${{PWD}}/{} {}".format(FILEDICT[fn], fn),
                            shell=True)
