using Base
using CSV
using XLSX
using DataFrames
using StatsBase
using Query
using Statistics

prefix = "/home/irina_pugach/aDNA_Indonesia_Guam/New_Guam/Data/for_Alex/Clean.HO-MOD.ancGuamALL.HO-ANC.YRI.FRE.Liangdao.MEGA"
# Read individual files
inds = DataFrame(CSV.File(string(prefix, ".ind");
                          header=["ind", "sex", "pop"], delim="\t"))
# Identify ancient samples
popinfo = DataFrame(CSV.File("/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/documentation/population_overview_SNPdata.csv";
                    header=true, delim="\t"))
ancient_pops = @from i in popinfo begin
               @where i.type == "a"
               @select i.pop
               @collect
           end

# Identify ancient samples in data set
ancient_indices = findall(x -> x in ancient_pops, inds["pop"])
ancient_individuals = inds[ancient_indices, "ind"]

# Read genotype file and only extract HGDP samples
geno_raw = map(collect, readlines(string(prefix, ".geno")))
geno = Array{Int}(undef, length(geno_raw), length(ancient_individuals))
Threads.@threads for i in 1:length(geno_raw[1])
    if inds[i, "ind"] in ancient_individuals
        k = findall(x -> x == inds[i, "ind"], ancient_individuals)[1]
        for j in 1:length(geno_raw)
            geno[j, k] = parse(Int, geno_raw[j][i])
        end
    end
end

# Pairwise differences
pwdiff = zeros(Float64, length(ancient_individuals), length(ancient_individuals))
numsites = zeros(Int64, length(ancient_individuals), length(ancient_individuals))
for i in 1:size(geno)[2]
    println(i)
    Threads.@threads for j in 1:size(geno)[2]
        if j > i
            println(join(["\t", string(j)]))
            notmissing = [a != 9 for a in geno[:,i]] .& [a != 9 for a in geno[:,j]]
            if sum(notmissing) > 0
                gts1 = [x == 1 ? sample([0, 2], 1) : x for x in geno[notmissing, i]]
                gts2 = [x == 1 ? sample([0, 2], 1) : x for x in geno[notmissing, j]]
                pwdiff[i, j] = sum(gts1 .!= gts2) / sum(notmissing)
                numsites[i, j] = sum(notmissing)
            end
        end
    end
end

# Save as DataFrame
pwdiff_df = DataFrame(pwdiff, convert.(String, ancient_individuals))
CSV.write("/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/analysis/qual/pairwise_differences/ancientSamples_pairwiseDiff.csv", pwdiff_df, delim = "\t")
numsites_df = DataFrame(numsites, convert.(String, ancient_individuals))
CSV.write("/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/analysis/qual/pairwise_differences/ancientSamples_numsites.csv", numsites_df, delim = "\t")

# Pairwise difference based on 1240K per window
geno_raw = map(collect, readlines("/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/analysis/genotypes/1240K/ancGuam.all.geno"))
geno = Array{Int}(undef, length(geno_raw), 2)
Threads.@threads for i in 1:length(geno_raw[1])
        for j in 1:length(geno_raw)
            geno[j, i] = parse(Int, geno_raw[j][i])
        end
end
guam_notmissing = [a != 9 for a in geno[:,1]] .& [a != 9 for a in geno[:,2]]
SP4210 = [x == 1 ? sample([0, 2], 1) : x for x in geno[guam_notmissing, 1]]
SP4211 = [x == 1 ? sample([0, 2], 1) : x for x in geno[guam_notmissing, 2]]
guam_1240K_pwdiff = sum(SP4210 .!= SP4211) / sum(guam_notmissing)
guam_1240K_nosites = sum(guam_notmissing)

snp_1240K = DataFrame(CSV.File("/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/analysis/genotypes/1240K/ancGuam.all.snp";
                               header=["snpname", "chr", "dummy", "pos", "ref", "alt"], delim="\t"))
notmissing_snps = snp_1240K[guam_notmissing, :]
notmissing_snps["prop"] = round.(Int, notmissing_snps["pos"] / 1e6)
numbersnps_per_window = by(notmissing_snps, [:chr, :prop], nrow)
notmissing_snps = join(notmissing_snps, numbersnps_per_window, on = [:chr, :prop])
notmissing_snps["window"] = string.(notmissing_snps["chr"], "_",  notmissing_snps["prop"])

guam_1240K_pwdiff_window = zeros(Float64, length(unique(notmissing_snps.window)))
guam_1240K_numsites_window = zeros(Int64, length(unique(notmissing_snps.window)))

windows = unique(notmissing_snps.window)
Threads.@threads for i in 1:length(windows)
    wd_indices = findall(x -> x .== windows[i], notmissing_snps.window)
    guam_1240K_pwdiff_window[i] = sum(SP4210[wd_indices] .!= SP4211[wd_indices]) / length(wd_indices)
    guam_1240K_numsites_window[i] = length(wd_indices)
end
guam_1240K_pwdiff_window_df = DataFrame(window = windows, diff = guam_1240K_pwdiff_window, nSites = guam_1240K_numsites_window)
guam_1240K_pwdiff_window_df["chr"] = [x[1] for x in split.(guam_1240K_pwdiff_window_df.window, "_")]

window_stats = by(guam_1240K_pwdiff_window_df[guam_1240K_pwdiff_window_df.nSites .>= 25, :], :chr, 
                  n = :diff => length,
                  mean = :diff => mean,
                  sd = :diff => std)
CSV.write("/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/results/QUAL_pairwiseDist_perwindow.csv",
          window_stats, delim="\t")
