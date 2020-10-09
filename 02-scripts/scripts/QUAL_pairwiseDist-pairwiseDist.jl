using Base
using CSV
using XLSX
using DataFrames
using StatsBase

# 1000Genomes on 1240K
## Read individual files
inds = DataFrame(CSV.File("../../03-data/reich_public_geno_v42.4.1240K.ind";
                          header=["ind", "sex", "pop"], delim="\t"))
## Information on 1000Genomes samples
sample_info = DataFrame(XLSX.readtable("../../01-documentation/1000Genomes_samplelist.xlsx", "Sample Info")...)
sample_names = convert.(String, sample_info["Sample"])

## Identify 1000Genomes samples in Reich data set
ind_names = [m != nothing ? m.captures[1] : nothing for m in match.(r"([A-Z]+[0-9]+).*", inds["ind"])]
tgenomes_samples = ind_names[in.(ind_names, [sample_names])]
tgenomes_indices = findall(x -> x in tgenomes_samples, ind_names)

## Read genotype file and only extract 1000Genomes samples
geno_raw = map(collect, readlines("/mnt/expressions/mateja/Early_modern_humans/David_files/v42.4.1240k/extracted.v42.4.1240K.geno"))
tgenomes_geno = Array{Int}(undef, length(geno_raw), length(tgenomes_indices))
Threads.@threads for i in 1:length(geno_raw[1])
    if ind_names[i] in tgenomes_samples
        k = findall(x -> x == ind_names[i], tgenomes_samples)[1]
        for j in 1:length(geno_raw)
            tgenomes_geno[j, k] = parse(Int, geno_raw[j][i])
        end
    end
end

## Pairwise differences
tgenomes_pwdiff = zeros(Float64, length(tgenomes_indices), length(tgenomes_indices))
for i in 1:size(tgenomes_geno)[2]
    println(i)
    Threads.@threads for j in 1:size(tgenomes_geno)[2]
        if j > i
            notmissing = [a != 9 for a in tgenomes_geno[:,i]] .& [a != 9 for a in tgenomes_geno[:,j]]
            if sum(notmissing) > 0
                gts1 = [x == 1 ? sample([0, 2], 1) : x for x in tgenomes_geno[notmissing, i]]
                gts2 = [x == 1 ? sample([0, 2], 1) : x for x in tgenomes_geno[notmissing, j]]
                tgenomes_pwdiff[i, j] = sum(gts1 .!= gts2) / sum(notmissing)
            end
        end
    end
end

## Save as DataFrame
tgenomes_pwdiff_df = DataFrame(tgenomes_pwdiff, convert.(String, tgenomes_samples))
CSV.write("/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/analysis/qual/pairwise_differences/1000Genomes_pairwisediff.csv", tgenomes_pwdiff_df, delim = "\t")

# HGDP on Human Origins
## Read individual files
inds = DataFrame(CSV.File("../../03-data/reich_public_geno_v42.4.1240K_HO.ind";
                          header=["ind", "sex", "pop"], delim="\t"))
## Information on 1000Genomes samples
sample_info = CSV.read("../../01-documentation/igsr-human_genome_diversity_project.tsv"; delim = "\t")
sample_names = convert.(String, sample_info["Sample name"])

## Identify 1000Genomes samples in Reich data set
ind_names = [m != nothing ? m.captures[1] : nothing for m in match.(r"([A-Z]+[0-9]+).*", inds["ind"])]
hgdp_samples = ind_names[in.(ind_names, [sample_names])]
hgdp_indices = findall(x -> x in hgdp_samples, ind_names)

## Read genotype file and only extract 1000Genomes samples
geno_raw = map(collect, readlines("../../03-data/extracted.v42.4.1240K_HO.geno"))
hgdp_geno = Array{Int}(undef, length(geno_raw), length(hgdp_indices))
Threads.@threads for i in 1:length(geno_raw[1])
    if ind_names[i] in hgdp_samples
        k = findall(x -> x == ind_names[i], hgdp_samples)[1]
        for j in 1:length(geno_raw)
            hgdp_geno[j, k] = parse(Int, geno_raw[j][i])
        end
    end
end

## Pairwise differences
hgdp_pwdiff = zeros(Float64, length(hgdp_indices), length(hgdp_indices))
for i in 1:size(hgdp_geno)[2]
    println(i)
    Threads.@threads for j in 1:size(hgdp_geno)[2]
        if j > i
            println(join(["\t", string(j)]))
            notmissing = [a != 9 for a in hgdp_geno[:,i]] .& [a != 9 for a in hgdp_geno[:,j]]
            if sum(notmissing) > 0
                gts1 = [x == 1 ? sample([0, 2], 1) : x for x in hgdp_geno[notmissing, i]]
                gts2 = [x == 1 ? sample([0, 2], 1) : x for x in hgdp_geno[notmissing, j]]
                hgdp_pwdiff[i, j] = sum(gts1 .!= gts2) / sum(notmissing)
            end
        end
    end
end

## Save as DataFrame
hgdp_pwdiff_df = DataFrame(hgdp_pwdiff[:,1:799], convert.(String, hgdp_samples[1:799]))
CSV.write("/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/analysis/qual/pairwise_differences/HGDP_pairwisediff.csv", hgdp_pwdiff_df, delim = "\t")
