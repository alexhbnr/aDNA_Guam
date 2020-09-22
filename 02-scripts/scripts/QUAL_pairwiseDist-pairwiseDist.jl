using Base
using CSV
using XLSX
using DataFrames
using StatsBase

# Read individual files
inds = DataFrame(CSV.File("../../03-data/reich_public_geno_v42.4.1240K.ind";
                          header=["ind", "sex", "pop"], delim="\t"))
# Information on HGDP samples
sample_info = DataFrame(XLSX.readtable("../../01-documentation/1000Genomes_samplelist.xlsx", "Sample Info")...)
sample_names = convert.(String, sample_info["Sample"])

# Identify HGDP samples in Reich data set
ind_names = [m != nothing ? m.captures[1] : nothing for m in match.(r"([A-Z]+[0-9]+).*", inds["ind"])]
hgdp_samples = ind_names[in.(ind_names, [sample_names])]
hgdp_indices = findall(x -> x in hgdp_samples, ind_names)

# Read genotype file and only extract HGDP samples
geno_raw = map(collect, readlines("/mnt/expressions/mateja/Early_modern_humans/David_files/v42.4.1240k/extracted.v42.4.1240K.geno"))
geno = Array{Int}(undef, length(geno_raw), length(hgdp_indices))
Threads.@threads for i in 1:length(geno_raw[1])
    if ind_names[i] in hgdp_samples
        k = findall(x -> x == ind_names[i], hgdp_samples)[1]
        for j in 1:length(geno_raw)
            geno[j, k] = parse(Int, geno_raw[j][i])
        end
    end
end

# Pairwise differences
pwdiff = zeros(Float64, length(hgdp_indices), length(hgdp_indices))
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
            end
        end
    end
end

# Save as DataFrame
pwdiff_df = DataFrame(pwdiff, convert.(String, hgdp_samples))
CSV.write("/mnt/genotyping/sk_pipelines/projects/aDNA_Guam/analysis/qual/pairwise_differences/HGDP_pairwisediff.csv", pwdiff_df, delim = "\t")
