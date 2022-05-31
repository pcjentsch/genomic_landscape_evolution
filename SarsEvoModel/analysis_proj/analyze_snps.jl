using VariantCallFormat
using PyCall
using Intervals
using DataFrames
using CSV
using Plots

const plotting_settings =
    (
        fontfamily="serif-roman",
        framestyle=:box,
        dpi=300,
        color_palette=palette(:seaborn_pastel)
    )
struct SNP
    ind::Int
    alt::String
    ref::Char
end
function import_binding_calc()
    py"""
    import sys
    sys.path.insert(0, "./SARS2_RBD_Ab_escape_maps")
    """
    bcalc = pyimport("bindingcalculator")
    b = bcalc.BindingCalculator(csv_or_url="SARS2_RBD_Ab_escape_maps/processed_data/escape_calculator_data.csv")
    return b
end
const b = import_binding_calc()
# https://www.nature.com/articles/s41598-021-99661-7#Sec1
const RBD = Interval(21560 + 331, 21560 + 531)


function get_snp(vcf::VariantCallFormat.Record)::Union{Nothing,SNP}
    ind = VCF.pos(vcf)
    ref = VCF.ref(vcf)
    alt = VCF.alt(vcf)
    if length(ref) == 1 && length(alt) == 1 && ref != "N"
        return SNP(
            ind,
            only(alt),
            only(ref),
        )
    else
        return nothing
    end
end

function Base.isequal(a::SNP, b::SNP)
    return a.ind == b.ind && a.ref == b.ref && a.alt == b.alt
end
function Base.hash(a::SNP)
    return hash((a.ind, a.ref, a.alt))
end
struct Genome
    snps::Set{SNP}
end
using ThreadsX
function parse_vcfs(vcf_path)
    genomes = ThreadsX.map(readdir(vcf_path)) do fname
        fpath = joinpath(vcf_path, fname)
        reader = VCF.Reader(open(fpath, "r"))
        records = collect(reader)
        meta = header(reader)
        close(reader)
        id = split(fname, ".") |> first
        snps = Iterators.filter(snp -> !isnothing(snp), Iterators.map(get_snp, records))
        genome = collect(snps)
        return (; id, genome)
    end |> DataFrame
    return genomes
end
using Dates
function get_data()
    genomes = parse_vcfs("vcf/")
    metadata_df = CSV.File("metadata_unzipped.tsv") |> DataFrame
    return genomes, metadata_df
end
using PyCall


function compute_binding_retained(genomes_w_metadata)
    inds_in_S_gene = filter(in(b.sites), map(snp -> snp.ind - 21560, g))
    return b.binding_retained(inds_in_S_gene)
end
function analyze(genomes, metadata)
    filter!(:date => d -> d ∉ ("2020", "2021"), metadata)
    metadata.id = replace.(metadata.strain, "/" => "_", "-" => "N")
    genomes_w_metadata = innerjoin(genomes, metadata; on=:id)
    genomes_w_metadata.date = Date.(genomes_w_metadata.date)
    genomes_w_metadata.binding_retained = compute_binding_retained(genomes_w_metadata)
    return genomes_w_metadata
end
function load_lineages()
    lineages = CSV.File("uk_lineages.csv") |> DataFrame
    lineages.id = replace.(lineages.taxon, "/" => "_", "-" => "N")
    return lineages
end

using Serialization
using StatsBase

function plot_samples(genomes_w_metadata)
    df_nonempty = filter(:genome => g -> !isempty(g), genomes_w_metadata)
    datelist = [count(==(d), df_nonempty.date) / count(==(d), genomes_w_metadata.date) for d in sort(unique(df_nonempty.date))]
    numlist = [count(==(d), df_nonempty.date) for d in sort(unique(df_nonempty.date))]
    p = plot(sort(unique(df_nonempty.date)), datelist; legend=false, xlabel="date", ylabel="fraction of samples", title="Fraction of UK samples with mutations in RBD", plotting_settings...)
    savefig(p, "rbd_fraction.png")
    p = plot(sort(unique(df_nonempty.date)), numlist; legend=false, xlabel="date", ylabel="no samples", title="number of UK samples by date", plotting_settings...)
    savefig(p, "num_samples.png")
    df = groupby(genomes_w_metadata, :date) |>
         df -> combine(df,
        :binding_retained => mean => :mean,
        :binding_retained => std => :std,
        :binding_retained => (x -> maximum(x)) => :upper_quartile,
        :binding_retained => (x -> minimum(x)) => :lower_quartile
    ) |> df -> sort(df, :date)

    p = plot(df.date, df.mean; ribbon=df.std)#(df.mean .- df.lower_quartile, df.upper_quartile .- df.mean))


    # return df

end

function dist(g1::Set{SNP}, g2::Set{SNP})
    d = 0
    for snp in g1
        if snp ∉ g2
            d += 1
        end
    end
    for snp in g2
        if snp ∉ g1
            d += 1
        end
    end
    return d#length(symdiff(g1.snps, g2.snps))
end
using MultivariateStats
function analyze_binding_retained(df)
    df.genome = Set.(df.genome)
    unique_genomes_df = unique(df, :genome)
    pairwise_distances = zeros(nrow(unique_genomes_df), nrow(unique_genomes_df))
    for i in 1:nrow(unique_genomes_df)
        binding_i = unique_genomes_df.binding_retained[i]#compute_binding_retained(unique_genomes_df.genome[i])
        genome_i::Set{SNP} = unique_genomes_df.genome[i]
        for j in 1:(i-1)
            binding_j = unique_genomes_df.binding_retained[j]#compute_binding_retained(unique_genomes_df.binding_retained[j])
            genome_j::Set{SNP} = unique_genomes_df.genome[j]
            snp_distance = dist(genome_i, genome_j)
            pairwise_distances[i, j] = snp_distance * ((1 - binding_i) + (1 - binding_j)) / 2
            pairwise_distances[j, i] = snp_distance * ((1 - binding_i) + (1 - binding_j)) / 2
        end
    end
    mds = fit(MDS, pairwise_distances; distances=true, maxoutdim=2)
    coord_matrix = predict(mds)
    unique_genomes_df.mds_x = coord_matrix[1, :]
    unique_genomes_df.mds_y = coord_matrix[2, :]
    return df, select(unique_genomes_df, [:genome, :mds_x, :mds_y])
end

function plot_mds(genomes_w_metadata, unique_genomes_df, lineages)
    df = innerjoin(genomes_w_metadata, unique_genomes_df; on=:genome) |> df ->
        innerjoin(lineages, df; on=:id)
    sort!(df, :date)
    p = scatter(unique_genomes_df.mds_x, unique_genomes_df.mds_y; label="samples", xlabel="MDS1", ylabel="MDS2",
        title="Approximate antigenic cartography w/ RBD, UK samples",
        markersize=1.5,
        markerstrokewidth=0.3,
        plotting_settings...)
    savefig(p, "multidimensional_scaling.png")

    anim = Animation()
    xlims = extrema(df.mds_x)
    ylims = extrema(df.mds_y)
    mds_x_accum = Float64[]
    mds_y_accum = Float64[]
    for (i, (key, gdf)) in enumerate(pairs(groupby(df, :date)))
        display(i)
        pts_df = unique(gdf, [:mds_x, :mds_y])
        append!(mds_x_accum, pts_df.mds_x)
        append!(mds_y_accum, pts_df.mds_y)
        p = scatter(mds_x_accum, mds_y_accum; label="samples", xlabel="MDS1", ylabel="MDS2",
            title="Estimated antigenic distance, $(key.date)",
            markersize=1.5,
            markerstrokewidth=0.3,
            xlims=xlims,
            ylims=ylims,
            plotting_settings...)
        frame(anim, p)
    end
    gif(anim, "animation.gif")
end
# function MDS(dataframe)
#     distance_matrix = 
# end


# genomes, metadata = get_data()
# genomes_w_metadata = analyze(genomes, metadata)
# serialize("genomes_metadata.dat", genomes_w_metadata)
# genomes_w_metadata = deserialize("genomes_metadata.dat")

# plot_samples(genomes_w_metadata)
# genomes_w_metadata, unique_genomes_df = analyze_binding_retained(genomes_w_metadata)
lineages_df = load_lineages()
plot_mds(genomes_w_metadata, unique_genomes_df, lineages_df)
# bcalc = pyimport("SARS2_RBD_Ab_escape_maps.bindingcalculator")
# b = bcalc.BindingCalculator(csv_or_url="processed_data/escape_calculator_data.csv")
