using PyCall
using Intervals
using DataFrames
using CSV
using StatsPlots
using ProgressMeter
using CodecZlib, TranscodingStreams, FASTX, BioSymbols, ThreadsX
using MultivariateStats
using Distributions, Random
using OnlineStats
using BenchmarkTools
using Infiltrator
using Cthulhu

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
    "Omicron (BA.5-like)" => "B.1.1.529.5", #omicron 4/5
    "Omicron (BA.4-like)" => "B.1.1.529.4", #omicron 4/5
    "Probable Omicron (BA.5-like)" => "B.1.1.529.5", #basal omicron
    "Probable Omicron (BA.4-like)" => "B.1.1.529.4", #basal omicron
    "Probable Omicron (XE-like)" => missing, 
    "Omicron (XE-like)" => missing, 
)
function parse_scorpio(scorpio)
    return scorpio_map[scorpio]
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
function make_snp_dict(recurrent_df, binding_sites)
    snps_inds = mapreduce(bsite -> findall(==(bsite), recurrent_df.ind), vcat, filter(in(binding_sites), recurrent_df.ind))
    top_n_snps = vcat(recurrent_df[1:100, :], recurrent_df[snps_inds, :]) |> unique
    snp_weight_dict = Dict(zip(top_n_snps.snp, top_n_snps.freq))
    return snp_weight_dict
end

function unique_genomes(df,snp_weight_dict)
    unique_genomes_df = copy(df)
    top_n_snps = Set(keys(snp_weight_dict))
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
    unique_genomes_df.filtered_genome = ThreadsX.map(genome -> filter(snp -> snp in top_n_snps, genome), unique_genomes_df.genome)
    dropmissing!(unique_genomes_df, [:closest_mapped_lineage, :mapped_lineage_position])
    return unique_genomes_df
end

function filter_genomes(unique_df)
    filtered_genome_df = deepcopy(unique_df)
    unique!(filtered_genome_df, [:closest_mapped_lineage, :filtered_genome])
    return filtered_genome_df
end
function measure_distance_multiplier(filter_df, snp_weight_dict)
    samples = nrow(filter_df)
    
    prog = Progress(samples)
    dm_weighted_snp = zeros(Float64, samples, samples)
    dm_mapped = zeros(Float64, samples, samples)
    function map_distances((i, j,))
       
        dm_weighted_snp[i, j] = d_snp
        dm_weighted_snp[j, i] = d_snp

        dm_mapped[i, j] = pt_d
        dm_mapped[j, i] = pt_d

        next!(prog)
    end
    ThreadsX.foreach(map_distances, lower_triangular(samples))
    # for (i,j) in lower_triangular(samples)
        
    # end

    return dm_weighted_snp,dm_mapped
end

function pairwise_distances(filtered_genome_df, snp_weight_dict)
    samples = nrow(filtered_genome_df)
    dm = zeros(Float64, samples, samples)
    prog = Progress(samples)
    snp_distance_weight = 5_000
    function map_distances((i, j,))
        genome_i = filtered_genome_df.filtered_genome[i]
        genome_j = filtered_genome_df.filtered_genome[j]
        mapped_lineage_position_i = filtered_genome_df.mapped_lineage_position[i]
        mapped_lineage_position_j = filtered_genome_df.mapped_lineage_position[j]
        binding_i = filtered_genome_df.binding_retained[i]
        binding_j = filtered_genome_df.binding_retained[j]

        snp_distance = weighted_dist(genome_i, genome_j, snp_weight_dict)
        pt_d = sqrt((mapped_lineage_position_i[1] - mapped_lineage_position_j[1])^2 + (mapped_lineage_position_i[2] - mapped_lineage_position_j[2])^2)
        binding_d = ((1 - binding_i) + (1 - binding_j)) / 2
        d = pt_d + snp_distance_weight * snp_distance * binding_d
        dm[i, j] = Float64(d)
        dm[j, i] = Float64(d)
        next!(prog)
    end
    ThreadsX.foreach(map_distances, lower_triangular(samples))

    return dm
end
function manifold_projection(dm, filter_genomes_df)
    n = nrow(filter_genomes_df)
    @assert all(size(dm) .== n)

    mds = fit(MDS, dm; distances=true, maxoutdim=2)
    coord_matrix = predict(mds)
    filter_genomes_df.mds_x = coord_matrix[1, :]
    filter_genomes_df.mds_y = coord_matrix[2, :]
    return mds
end
function pairwise_distances_bootstrap_ci(unique_genomes_df, filtered_genome_df, snp_weight_dict, map_distances_fn)
    samples = nrow(filtered_genome_df)
    dm = zeros(Float32, samples, samples)
    prog = Progress(samples)
    function map_distances((i, j,))
        d = map_distances_fn(
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
    noise_dist = LogNormal(0.0, 0.25)

    filtered_genome_df.mds_x_ci = [Variance() for _ in 1:samples]
    filtered_genome_df.mds_y_ci = [Variance() for _ in 1:samples]
    noise_vec = zeros(div(samples * (samples - 1), 2))
    no_bootstraps = 10
    for k in 1:no_bootstraps
        rand!(noise_dist, noise_vec)
        ind = 1
        for i in 1:samples
            for j in 1:i
                dm[i, j] += noise_vec[ind]
                dm[i, j] += noise_vec[ind]
            end
        end
        display(dm)
        mds = fit(MDS, dm; distances=true, maxoutdim=2)
        coord_matrix = predict(mds)
        fit!.(filtered_genome_df.mds_x_ci, coord_matrix[1, :])
        fit!.(filtered_genome_df.mds_y_ci, coord_matrix[2, :])
    end

    for r in eachrow(filtered_genome_df)
        inds = findall(==(r.filtered_genome), unique_genomes_df.filtered_genome)
        unique_genomes_df.mds_x[inds] .= r.mds_x
        unique_genomes_df.mds_y[inds] .= r.mds_y
    end
    return unique_genomes_df, filtered_genome_df, dm
end


