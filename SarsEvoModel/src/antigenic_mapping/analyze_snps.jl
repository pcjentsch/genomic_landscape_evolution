using PyCall
using Intervals
using DataFrames
using CSV
using StatsPlots
using ProgressMeter
using CodecZlib, TranscodingStreams, FASTX, BioSymbols, ThreadsX
using MultivariateStats

const gene_index = [
    Interval(266, 21555),
    Interval(21563, 25384),
    Interval(28274, 29533),
    Interval(25393, 26220),
    Interval(26523, 27191),
    Interval(27894, 28259),
    Interval(27394, 27759),
    Interval(27756, 27887),
    Interval(27202, 27387),
    Interval(26245, 26472),
    Interval(29558, 29674),
]
const gene_names = [
    "orf1ab",
    "S",
    "N",
    "ORF3a",
    "M",
    "ORF8",
    "ORF7a",
    "ORF7b",
    "ORF6",
    "E",
    "ORF10",
]
const py_bcalc = PyNULL()
# https://www.nature.com/articles/s41598-021-99661-7#Sec1
const S_gene_ind = 21560
const RBD = Interval(S_gene_ind + 331, S_gene_ind + 531)


function __init__()
    path = joinpath(@__DIR__, "../../deps/SARS2_RBD_Ab_escape_maps/")
    display(abspath(path))
    py"""
    import sys
    sys.path.insert(0,$path)
    """
    bcalc = pyimport("bindingcalculator")
    copy!(py_bcalc, bcalc.BindingCalculator(csv_or_url=joinpath(@__DIR__, "../../deps/SARS2_RBD_Ab_escape_maps/processed_data/escape_calculator_data.csv")))

end

function compute_binding_retained(g)
    return py_bcalc.binding_retained(
        [(snp.ind - S_gene_ind) for snp in g])
end

function analyze_binding_retained(df, recurrent_df)
    snps_inds = mapreduce(bsite -> findall(==(bsite), recurrent_df.ind), vcat, filter(in(binding_sites), recurrent_df.ind))
    top_n_snps = vcat(recurrent_df[1:100, :], recurrent_df[snps_inds, :]) |> unique
    top_n_snps_set = Set(top_n_snps.snp)
    snp_weight_dict = Dict(zip(top_n_snps.snp, top_n_snps.freq))
    df.filtered_genome = ThreadsX.map(genome -> Set(Iterators.filter(snp -> snp in top_n_snps_set, genome)), df.genome)
    unique_genomes_df = unique(df, :filtered_genome)
    samples = nrow(unique_genomes_df)
    @info samples
    unique_genomes_subsampled = unique_genomes_df[sample(1:nrow(unique_genomes_df), samples; replace=false), :]
    pairwise_distances = zeros(samples, samples)
    for i in 1:samples
        binding_i = unique_genomes_df.binding_retained[i]
        genome_i::Set{SNP} = unique_genomes_df.filtered_genome[i]
        for j in 1:(i-1)
            binding_j = unique_genomes_df.binding_retained[j]
            genome_j::Set{SNP} = unique_genomes_df.filtered_genome[j]
            snp_distance = dist(genome_i, genome_j, snp_weight_dict)
            pairwise_distances[i, j] = snp_distance * ((1 - binding_i) + (1 - binding_j)) / 2
            pairwise_distances[j, i] = snp_distance * ((1 - binding_i) + (1 - binding_j)) / 2
        end
    end
    mds = fit(MDS, pairwise_distances; distances=true, maxoutdim=2)
    coord_matrix = predict(mds)
    unique_genomes_subsampled.mds_x = coord_matrix[1, :]
    unique_genomes_subsampled.mds_y = coord_matrix[2, :]
    # unique_genomes_df.mds_z = coord_matrix[3, :]
    return df, unique_genomes_subsampled, mds
end



# gene_occurrence = groupby(recurrent_df,:gene_name) |> df -> combine(df, nrow) |> df -> filter(:gene_name => n ->length(n)==1,df)
# gene_occurrence.freq = gene_occurrence.occurrence_sum .- map(i -> abs(first(i) - last(i)),gene_index)
# #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7161481/
# const gene_lengths = [




#     for (key,gdf) in pairs(groupby(recurrent_df,:gene_name))
#         savefig(plot(gdf.occurrence; title = key.gene_name),"$(key.gene_name).png")
#     end


# gcf = CSV.File("GCF_009858895.2_ASM985889v3_genomic.gff"; delim="\t") |> DataFrame |> 
#     df-> dropmissing(df,:Column9)
# parse_dict(s) = Dict(map(t -> Pair(t...), split.(split(s, ";"), "=")))
# gcf.keys = parse_dict.(gcf.Column9)
# parse_line(l) = occursin("mature", l.Column3) ? l.keys["ID"] : l.keys["product"]
# function parse_line(l) 
#     if haskey(l.keys, "product") 
#         l.keys["ID"] 
#     elseif haskey(l.keys,"ID")

# top 10 most mutated snps

# 10×6 DataFrame
#  Row │ ID       occurrence  snp               freq    ⋯
#      │ String7  Int64       SNP               Float64 ⋯
# ─────┼─────────────────────────────────────────────────
#    1 │ G7328T         6766  SNP(0x1ca0, 'T')  0.00272 ⋯
#    2 │ T27752C        6583  SNP(0x6c68, 'C')  0.00265
#    3 │ T7124C         5086  SNP(0x1bd4, 'C')  0.00205
#    4 │ C27638T        4746  SNP(0x6bf6, 'T')  0.00191
#    5 │ T29402G        4744  SNP(0x72da, 'G')  0.00191 ⋯
#    6 │ T3037C         4206  SNP(0x0bdd, 'C')  0.00169
#    7 │ T28095A        3584  SNP(0x6dbf, 'A')  0.00144
#    8 │ A10323G        3255  SNP(0x2853, 'G')  0.00131
#    9 │ A21137G        3130  SNP(0x5291, 'G')  0.00126 ⋯
#   10 │ G21618C        3006  SNP(0x5472, 'C')  0.00121

#T3037C has some results 
#e.g. https://journals.plos.org/Plospathogens/article?id=10.1371/journal.ppat.1009849
# and https://www.nature.com/articles/s41586-021-03610-3 recognize it as a common SNP and iSNV

# A10323G is mentioned as common in
# https://web.archive.org/web/20210204225512id_/https://www.medrxiv.org/content/medrxiv/early/2021/01/09/2021.01.08.21249379.full.pdf
# smooth out orf3a plot
# look at k-means of distance matrix