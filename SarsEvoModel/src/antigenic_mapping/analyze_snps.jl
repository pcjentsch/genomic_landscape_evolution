using PyCall
using Intervals
using DataFrames
using CSV
using StatsPlots
using ProgressMeter
using CodecZlib, TranscodingStreams, FASTX, BioSymbols, ThreadsX
using MultivariateStats
const scorpio_map = Dict(
    missing => missing,
    "Alpha (B.1.1.7-like)" => "B.1.1.7",
    "Delta (B.1.617.2-like)" => "B.1.617.2",
    "Epsilon (B.1.429-like)" => "B.1.429",
    "Epsilon (B.1.427-like)" => "B.1.429",
    "Zeta (P.2-like)" => "B.1.1.28.1", #like gamma
    "Iota (B.1.526-like)" => "B.1.526+E484K", #the closest we have to iota
    "Omicron (BA.1-like)" => "B.1.1.529",
    "Probable Omicron (BA.1-like)" => "B.1.1.529",
    "Delta (AY.4-like)" => "B.1.617.2", #basal delta I guess
    "Gamma (P.1-like)" => "B.1.1.28.1",
    "Beta (B.1.351-like)" => "B.1.351",
    "Eta (B.1.525-like)" => missing,
    "B.1.1.318-like" => missing,
    "A.23.1-like" => missing,
    "Theta (P.3-like)" => "B.1.1.28.1", #like gamma
    "Lambda (C.37-like)" => "B.1.1.1.37",
    "B.1.1.7-like+E484K" => "B.1.1.7+E484K",
    "B.1.617.1-like" => "B.1.617.1",
    "B.1.617.3-like" => "B.1.617.1",
    "Mu (B.1.621-like)" => "B.1.621",
    "Delta (B.1.617.2-like) +K417N" => "B.1.617.2(AY.2)+K417N",
    "Delta (AY.4.2-like)" => "B.1.617.2", #basal delta I guess
    "Omicron (BA.2-like)" => "B.1.1.529", #basal omicron
    "Omicron (Unassigned)" => "B.1.1.529", #basal omicron
    "Probable Omicron (Unassigned)" => "B.1.1.529", #basal omicron
    "Probable Omicron (BA.2-like)" => "B.1.1.529", #basal omicron
    "Omicron (BA.3-like)" => "B.1.1.529", #basal omicron
    "Probable Omicron (BA.3-like)" => "B.1.1.529", #basal omicron
    "Omicron (BA.5-like)" => "B.1.1.529", #basal omicron
    "Omicron (BA.4-like)" => "B.1.1.529", #basal omicron
    "Probable Omicron (BA.5-like)" => "B.1.1.529", #basal omicron
    "Probable Omicron (BA.4-like)" => "B.1.1.529", #basal omicron
    "Probable Omicron (XE-like)" => "B.1.1.529", #basal omicron
    "Omicron (XE-like)" => missing, #basal omicron
)
function parse_scorpio(scorpio)
    return scorpio_map[scorpio]
    # "B.1.617.2(AY.1)+K417N"
    # "B.1.617.2(AY.3)+E484Q"
    # "D614G"
end


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

function unique_genomes(df, recurrent_df, binding_sites)
    snps_inds = mapreduce(bsite -> findall(==(bsite), recurrent_df.ind), vcat, filter(in(binding_sites), recurrent_df.ind))
    top_n_snps = vcat(recurrent_df[1:150, :], recurrent_df[snps_inds, :]) |> unique
    snp_weight_dict = Dict(zip(top_n_snps.snp, top_n_snps.freq))
    unique_genomes_df = deepcopy(df)
    unique!(unique_genomes_df, :genome)
    unique_genomes_df.closest_mapped_lineage = Vector{Union{String,Missing}}(undef, nrow(unique_genomes_df))
    unique_genomes_df.mapped_lineage_position = Vector{Union{Tuple{Float64,Float64},Missing}}(undef, nrow(unique_genomes_df))
    for r in eachrow(unique_genomes_df)
        map_key = parse_scorpio(r.scorpio_call)
        map_position = antigenic_map_paper[map_key]
        r.closest_mapped_lineage = map_key
        r.mapped_lineage_position = ismissing(map_position) ? missing : map_position[1:2]
    end
    @show count(ismissing, unique_genomes_df.closest_mapped_lineage) / nrow(unique_genomes_df)
    dropmissing!(unique_genomes_df, [:closest_mapped_lineage, :mapped_lineage_position])

    unique_genomes_df.filtered_genome = ThreadsX.map(genome -> filter(snp -> snp in top_n_snps.snp, genome), unique_genomes_df.genome)
    filtered_genome_df = deepcopy(unique_genomes_df)
    unique!(filtered_genome_df, :filtered_genome)
    return unique_genomes_df, filtered_genome_df, snp_weight_dict
end
using BenchmarkTools
using Infiltrator
using Cthulhu
function compute_antigenic_distance(
    mapped_lineage_position_i::Tuple{Float64,Float64},
    mapped_lineage_position_j::Tuple{Float64,Float64},
    genome_i,
    genome_j,
    binding_i::Float64,
    binding_j::Float64,
    snp_weight_dict
)
    snp_distance = weighted_dist(genome_i, genome_j, snp_weight_dict)
    pt_d = (mapped_lineage_position_i[1] - mapped_lineage_position_j[1])^2 + (mapped_lineage_position_i[2] - mapped_lineage_position_j[2])^2
    d = pt_d + ((1 - binding_i) + (1 - binding_j)) / 2 #snp_distance * 1000#
    return d
end
function pairwise_distances(unique_genomes_df, filtered_genome_df, snp_weight_dict)
    samples = nrow(filtered_genome_df)
    dm = zeros(Float32, samples, samples)
    prog = Progress(samples)
    function map_distances((i, j,))
        d = compute_antigenic_distance(
            filtered_genome_df.mapped_lineage_position[i],
            filtered_genome_df.mapped_lineage_position[j],
            filtered_genome_df.filtered_genome[i],
            filtered_genome_df.filtered_genome[j],
            filtered_genome_df.binding_retained[i],
            filtered_genome_df.binding_retained[j],
            snp_weight_dict,
        )
        dm[i, j] = Float32(d)
        dm[j, i] = Float32(d)
        next!(prog)
    end
    ThreadsX.foreach(map_distances, lower_triangular(samples))
    mds = manifold_projection(dm, filtered_genome_df)
    select!(filtered_genome_df, [:genome, :filtered_genome, :mds_x, :mds_y])
    unique_genomes_df.mds_x = zeros(nrow(unique_genomes_df))
    unique_genomes_df.mds_y = zeros(nrow(unique_genomes_df))
    for r in eachrow(filtered_genome_df)
        inds = findall(==(r.filtered_genome), unique_genomes_df.filtered_genome)
        unique_genomes_df.mds_x[inds] .= r.mds_x
        unique_genomes_df.mds_y[inds] .= r.mds_y
    end
    return unique_genomes_df, filtered_genome_df, dm
end
function manifold_projection(dm, unique_genomes_df)
    n = nrow(unique_genomes_df)
    @assert all(size(dm) .== n)

    mds = fit(MDS, dm; distances=true, maxoutdim=2)
    coord_matrix = predict(mds)
    unique_genomes_df.mds_x = coord_matrix[1, :]
    unique_genomes_df.mds_y = coord_matrix[2, :]
    return mds
end

# expectation maximization

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