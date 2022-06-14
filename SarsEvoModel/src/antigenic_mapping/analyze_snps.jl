using VariantCallFormat
using PyCall
using Intervals
using DataFrames
using CSV
using Plots
using ProgressMeter
using CodecZlib, TranscodingStreams, FASTX, BioSymbols, ThreadsX

const plotting_settings =
    (
        fontfamily="serif-roman",
        framestyle=:box,
        dpi=300,
        color_palette=palette(:seaborn_pastel)
    )



function SNPs_from_fastas(aligned_fastas_dir::String)
    reference_reader = FASTA.Reader(open("../reference.fasta"))
    reference = FASTA.sequence(only(reference_reader))
    files = readdir(aligned_fastas_dir) |>
            l -> filter(s -> occursin("aligned", s), l) |>
                 l -> map(s -> joinpath(aligned_fastas_dir, s), l)
    display(files)
    file_streams = map(FASTA.Reader ∘ GzipDecompressorStream ∘ open, files)
    undef_base = gap(DNA)
    SNPs_by_sample = ThreadsX.map(Iterators.flatten(file_streams)) do sample_record
        id = FASTA.identifier(sample_record)
        genome = SNP[]
        sample_sequence = FASTA.sequence(sample_record)
        for (ind, (sample_base, ref_base)) in enumerate(zip(sample_sequence, reference))
            if ref_base != undef_base && sample_base != undef_base && ref_base != sample_base
                snp = SNP(ref_base, ind, sample_base)
                push!(genome, snp)
            end
        end
        return (; id, genome)
    end
    foreach(close, file_streams)
    return SNPs_by_sample
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


#matutils in usher tree
#number of independent acquisitions 

const b = import_binding_calc()
# https://www.nature.com/articles/s41598-021-99661-7#Sec1
const S_gene_ind = 21560
const binding_sites = Set(py"""list($(b.sites))""") .+ S_gene_ind
const RBD = Interval(S_gene_ind + 331, S_gene_ind + 531)


function get_snp(vcf::VariantCallFormat.Record)::Union{Nothing,SNP}
    ind = VCF.pos(vcf)
    ref = VCF.ref(vcf)
    alt = VCF.alt(vcf)

    if length(ref) == 1 && length(alt) == 1
        return SNP(
            only(ref),
            ind,
            only(only(alt)),
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

function parse_vcfs(vcf_path)
    vcf_files = readdir(vcf_path)
    p = Progress(length(vcf_files), 1)
    genomes = ThreadsX.map(readdir(vcf_path)) do fname
        next!(p)
        fpath = joinpath(vcf_path, fname)
        reader = VCF.Reader(open(fpath, "r"))
        records = collect(reader)
        meta = header(reader)
        close(reader)
        id = split(fname, ".") |> first
        snps = Iterators.filter(snp -> !isnothing(snp), Iterators.map(get_snp, records))
        genome = Set(collect(snps))
        return (; id, genome)
    end |> DataFrame
    return genomes
end

using Dates
function get_data_vcf()
    genomes = parse_vcfs("vcf/")
    metadata_df = CSV.File("metadata_unzipped.tsv") |> DataFrame
    filter!(:date => d -> d ∉ ("2020", "2021"), metadata)
    metadata.id = replace.(metadata.strain, "/" => "_", "-" => "N")
    genomes_w_metadata = innerjoin(genomes, metadata; on=:id)
    genomes_w_metadata.date = Date.(genomes_w_metadata.date)
    genomes_w_metadata.in_rbd = map(genomes_w_metadata.genome) do snps
        return Set(snp for snp in snps if snp.ind in binding_sites)
    end
    unique_genomes = select(unique(genomes_w_metadata, :in_rbd), [:in_rbd])
    unique_genomes.binding_retained = @showprogress map(compute_binding_retained, unique_genomes.in_rbd)
    genomes_w_metadata = innerjoin(unique_genomes, genomes_w_metadata; on=:in_rbd)
    return genomes_w_metadata
end
function get_data_fasta(fasta_path, metadata_path)
    snps = SNPs_from_fastas(fasta_path) |> DataFrame
    metadata = CSV.File(metadata_path) |> DataFrame
    filter!(:Collection_Date => d -> d ∉ ("2020", "2021"), metadata)

    genomes_w_metadata = innerjoin(snps, metadata; on=:id => :strain)
    genomes_w_metadata.date = Date.(genomes_w_metadata.Collection_Date)
    genomes_w_metadata.in_rbd = map(genomes_w_metadata.genome) do snps
        return Set(snp for snp in snps if snp.ind in binding_sites)
    end
    unique_genomes = select(unique(genomes_w_metadata, :in_rbd), [:in_rbd])
    unique_genomes.binding_retained = @showprogress map(compute_binding_retained, unique_genomes.in_rbd)
    genomes_w_metadata = innerjoin(unique_genomes, genomes_w_metadata; on=:in_rbd)
    return genomes_w_metadata
end
using PyCall


function compute_binding_retained(g)
    return b.binding_retained(
        [(snp.ind - S_gene_ind) for snp in g])
end

function analyze(genomes, metadata)

end

function load_lineages()
    lineages = CSV.File("uk_lineages.csv") |> DataFrame
    lineages.id = replace.(lineages.taxon, "/" => "_", "-" => "N")
    lineages.scorpio_call = map(lineages.scorpio_call) do scorpio_call
        ismissing(scorpio_call) && return missing
        if occursin("Delta", scorpio_call)
            return "delta"
        else
            return "everything else"
        end
    end
    return lineages
end


function plot_samples(genomes_w_metadata)
    df_nonempty = filter(:in_rbd => g -> !isempty(g), genomes_w_metadata)
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

    p = plot(df.date, df.mean; ribbon=df.std)
end

function dist(g1::Set{SNP}, g2::Set{SNP}, SNP_weights)
    d = 0
    for snp in g1
        if snp ∉ g2
            d += SNP_weights[snp]
        end
    end
    for snp in g2
        if snp ∉ g1
            d += SNP_weights[snp]
        end
    end
    return d#length(symdiff(g1.snps, g2.snps))
end
using MultivariateStats

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
# sooo the outliers appear to be the ones with a low binding..


function get_snp(line::T) where {T<:AbstractString}
    ref = line[1]
    alt = line[end]
    ind = parse(Int, line[2:end-1])
    return SNP(ref, ind, alt)
end

function parse_recurrent_mutations()
    recurrent = CSV.File("filtered_mutations.tsv") |> DataFrame
    recurrent.snp = map(get_snp, recurrent.ID)
    recurrent.freq = recurrent.occurrence ./ sum(recurrent.occurrence)
    sort!(recurrent, :freq; rev=true)
    name_by_ind(ind) = gene_names[findall(interval -> ind in interval, gene_index)]
    recurrent.ind = map(snp -> snp.ind, recurrent.snp)
    recurrent.gene_name = map(name_by_ind, recurrent.ind)
    filter!(:gene_name => n -> length(n) == 1, recurrent)
    recurrent.gene_name = only.(recurrent.gene_name)
    return recurrent
end

function plot_mds(genomes_df, unique_df, lineages)
    # df = innerjoin(genomes_df, unique_df; on=:filtered_genome => :filtered_genome) |> df ->
    #     innerjoin(lineages, df; on=:id)
    # sort!(df, :date)
    unique_df = innerjoin(unique_df, lineages; on=:id => :taxon, makeunique=true)

    p = plot()
    display(unique_df)
    # unique_df.is_omicron = map(s -> occursin("BA", s), unique_df.Pangolin)
    for (key, gdf) in pairs(groupby(unique_df, :scorpio_call))
        p = scatter!(p, gdf.mds_x, gdf.mds_y; label=key.scorpio_call, xlabel="MDS1", ylabel="MDS2",
            title="Approximate antigenic cartography w/ RBD, UK samples",
            markersize=1.5,
            markerstrokewidth=0.3,
            size=(800, 500),
            plotting_settings...)
    end
    savefig(p, "multidimensional_scaling.png")

    # anim = Animation()
    # xlims = extrema(df.mds_x)
    # ylims = extrema(df.mds_y)
    # mds_x_accum = Float64[]
    # mds_y_accum = Float64[]
    # for (i, (key, gdf)) in enumerate(pairs(groupby(df, :date)))
    #     display(i)
    #     pts_df = unique(gdf, [:mds_x, :mds_y])
    #     append!(mds_x_accum, pts_df.mds_x)
    #     append!(mds_y_accum, pts_df.mds_y)
    #     p = scatter(mds_x_accum, mds_y_accum; label="samples", xlabel="MDS1", ylabel="MDS2",
    #         title="Estimated antigenic distance, $(key.date)",
    #         markersize=1.5,
    #         markerstrokewidth=0.3,
    #         xlims=xlims,
    #         ylims=ylims,
    #         plotting_settings...)
    #     frame(anim, p)
    # end
    # gif(anim, "animation.gif")
end

datasets = (
    ("uk", "2021_11_28_gisaid_ncov_ingest/alignments/", "2021_11_28_gisaid_ncov_ingest/uk_sequences_metadata_new.tsv"),
    ("usa", "../newdata/alignments", "../newdata/usa_sequences_metadata.tsv")
)
# name, alignments, metadata = datasets[1]
# genome_set = map([datasets[2]]) do ()
# genomes_w_metadata = get_data_fasta(alignments, metadata)
#  serialize("genomes_metadata.dat", genomes_w_metadata)

# genomes_w_metadata = deserialize("genomes_metadata.dat")
recurrent_df = parse_recurrent_mutations()

plot_samples(genomes_w_metadata)
genomes_w_metadata_new, unique_genomes_df_new, mds = analyze_binding_retained(genomes_w_metadata, recurrent_df)
lineages_df = load_lineages()
plot_mds(genomes_w_metadata_new, unique_genomes_df_new, lineages_df)
# end

function orfplot()
    orf3a_coords = [
        ("Topological domain, Extracellular", Interval(1, 126)),
        ("Transmembrane, Helical", Interval(127, 183)),
        ("Topological domain, Cytoplasmic", Interval(184, 201)),
        ("Transmembrane, Helical", Interval(202, 279)),
        ("Topological domain, Extracellular", Interval(280, 303)),
        ("Transmembrane, Helical", Interval(304, 378)),
        ("Topological domain, Cytoplasmic", Interval(379, 825)),
    ]

    orf3a_df = groupby(filter(:gene_name => ==("ORF3a"), recurrent_df), :ind) |> df -> combine(df, :occurrence => sum)
    orf3a_df.description = map(eachrow(orf3a_df)) do r
        inds = filter(i -> (r.ind - 25393) in last(i), orf3a_coords)
        isempty(inds) && return missing
        return inds |> only |> first
    end
    dropmissing!(orf3a_df)
    orf3a_plot = plot()
    for (key, gdf) in pairs(groupby(orf3a_df, :description))
        scatter!(orf3a_plot, gdf.ind, gdf.occurrence_sum; label=key.description,
            xlabel="nucleotide index", ylabel="total recurrences",
            title="orf3a recurrences by index", markersize=2, markerstrokewidth=0.1, plotting_settings...)
    end
    savefig(orf3a_plot, "orf3a_scatter_2.png")
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