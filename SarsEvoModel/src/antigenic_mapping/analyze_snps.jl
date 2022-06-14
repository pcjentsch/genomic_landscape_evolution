using PyCall
using Intervals
using DataFrames
using CSV
using Plots
using ProgressMeter
using CodecZlib, TranscodingStreams, FASTX, BioSymbols, ThreadsX
using MultivariateStats

const plotting_settings =
    (
        fontfamily="serif-roman",
        framestyle=:box,
        dpi=300,
        color_palette=palette(:seaborn_pastel)
    )


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


#matutils in usher tree
#number of independent acquisitions 

const b = import_binding_calc()
# https://www.nature.com/articles/s41598-021-99661-7#Sec1
const S_gene_ind = 21560
const binding_sites = Set(py"""list($(b.sites))""") .+ S_gene_ind
const RBD = Interval(S_gene_ind + 331, S_gene_ind + 531)



function compute_binding_retained(g)
    return b.binding_retained(
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
# smooth out orf3a plot
# look at k-means of distance matrix