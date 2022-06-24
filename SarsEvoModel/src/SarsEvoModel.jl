module SarsEvoModel
export main, compute_antigenic_mapping_from_samples
using BenchmarkTools

#Ontario population
const initial_pop = 67.22e6
#Size of antigenic grid
const w = 25
const h = 25
function map_coords_to_model_space(coords_x, coords_y)
    return trunc.(Int, (coords_x * (w / 10), coords_y * (h / 10)))
end
include("antigenic_mapping/antigenic_map.jl")
include("antigenic_mapping/genomic_types.jl")

include("model.jl")
include("data.jl")
include("plotting.jl")
include("antigenic_mapping/analyze_snps.jl")
include("antigenic_mapping/antigenic_map_plotting.jl")
include("antigenic_mapping/io.jl")

lower_triangular(n) = ((i, j) for i in 1:n, j in 1:n if j < i)

function main()

    third_wave_begin = Date(2021, 03, 01)
    third_wave_end = Date(2021, 08, 01)
    location_data = UKLocationData()

    inv_infectious_period = 1 / 10
    r_0 = 2.5
    transmission_rate = r_0 * inv_infectious_period

    β = Float64[transmission_rate for i in 1:w, j in 1:h]
    sigma_matrix = Float64[sigma(i - k, j - l) for i in 1:w, j in 1:h, k in 1:w, l in 1:h]

    params = ModelParameters(
        β,
        1 / 10,
        0.07,
        initial_pop,
        0.00001,
        sigma_matrix,
    )
    sol = run(location_data, third_wave_begin, third_wave_end, params)
    plot_solution(sol)

    # # plot_antigenic_map()
    # return location_data
end

using Clustering, KernelDensity, AverageShiftedHistograms
function make_antigenic_map()
    binding_sites = Set(py"""list($(py_bcalc.sites))""") .+ S_gene_ind
    datasets = (
        (
            "uk",
            datapath("uk_sequences/alignments/"),
            datapath("uk_sequences/uk_sequences_metadata_new.tsv"),
            datapath("uk_sequence/uk_lineages.csv"),
        ),
        (
            "usa",
            datapath("usa_sequences/alignments/"),
            datapath("usa_sequences/usa_sequences_metadata.tsv"),
            datapath("usa_sequences/usa_lineages.csv"),
        )
    )
    for dataset in datasets[2:2]
        name, alignments, metadata, lineage_path = dataset
        recurrent_df = parse_recurrent_mutations(datapath("filtered_mutations.tsv"))
        genomes_w_metadata = serial_load(
            () -> get_data_fasta(alignments, metadata, binding_sites, lineage_path),
            datapath("genomes_$name.data")
        )
        unique_df, snp_weight_dict = unique_genomes(genomes_w_metadata, recurrent_df, binding_sites)
        unique_df, filter_df = pairwise_distances(unique_df, snp_weight_dict)
        plot_mds(name, unique_df)

        # p = scatter(unique_df.mds_x, unique_df.mds_y)
        # for gdf in groupby(unique_df,:)
        unique_df.week_submitted = round.(unique_df.date_submitted, Week)

        bounds_x = extrema(unique_df.mds_x)
        bounds_y = extrema(unique_df.mds_y)
        x_grid = LinRange(bounds_x..., 1000)
        y_grid = LinRange(bounds_y..., 1000)
        anim = Animation()
        grid = zeros(1000, 1000)
        sort!(unique_df, :week_submitted)
        for (key, gdf) in pairs(groupby(unique_df, :week_submitted))
            # kde_data = kde((gdf.mds_x, gdf.mds_y); npoints=(128, 128))
            o = ash(gdf.mds_x, gdf.mds_y; rngx=x_grid, rngy=y_grid, mx=100, my=100)
            # for ind in CartesianIndices(grid)
            #     i, j = Tuple(ind)
            #     grid[ind] = pdf(o, x_grid[i], y_grid[j])
            # end
            p = plot()
            scatter!(p, gdf.mds_x, gdf.mds_y; xlims=bounds_x, ylims=bounds_y, markersize=1.5,
                markerstrokewidth=0.3,
                size=(400, 300),
                plotting_settings...)
            plot!(p, o.rngx, o.rngy, o.z; title=key.date_submitted, xlims=bounds_x, ylims=bounds_y, plotting_settings...)

            # p = heatmap(o; plotting_settings...)
            frame(anim, p)
        end
        gif(anim, plots_path("$(name)_mds_density"; filetype="gif"))
        return unique_df, filter_df
    end
end

# f = Optimization.OptimizationFunction(clust_objective, Optimization.AutoForwardDiff())
# prob = Optimization.OptimizationProblem(clust_objective, initial_weight_vec, p;
#     lb=[0.0 for _ in initial_weight_vec], ub=[1.0 for _ in initial_weight_vec])
# sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited())
# return sol
# return unique_df, dm
function clust_objective(snp_weight_vec, p)
    snp_dict, unique_df = p
    snp_dict = Dict(zip(keys(snp_dict), snp_weight_vec))
    dm = pairwise_distances(unique_df, snp_dict)
    clust = hclust(dm; linkage=:single)
    nclusters = 3:8
    perf = map(nclusters) do k
        unique_df.clusters = cutree(clust; k)
        cat_df = groupby(unique_df, [:clusters, :category]) |>
                 gdf -> combine(gdf, nrow) |> df -> sort(df, :clusters)
        display(cat_df)

    end
    unique_df.clusters = cutree(clust; k=nclusters[argmax(perf)])
    cat_df = groupby(unique_df, [:clusters, :category]) |>
             gdf -> combine(gdf, nrow) |> df -> sort(df, :clusters)
    display(cat_df)
    return Clustering.randindex(cutree(clust; k=4), unique_df.category)[1]# maximum(perf)
end

function analyze_recurrences()
    recurrent_df = parse_recurrent_mutations(joinpath(@__DIR__, "../data/filtered_mutations.tsv"))
    reference_reader = FASTA.Reader(open(joinpath(@__DIR__, "../data/reference.fasta")))
    reference = FASTA.sequence(only(reference_reader))[1:end-2]
    n = 0
    df1 = (filter([:snp, :ref] => (s, r) -> reference[s.ind] != convert(DNA, r), recurrent_df))
    recurrent_df.actual_ref = map(s -> reference[s.ind], recurrent_df.snp)
    df2 = filter([:ref, :actual_ref] => (r, ra) -> ra != convert(DNA, r), recurrent_df)
    return df1, df2
    # for r in eachrow(recurrent_df)
    #     (; ind, alt) = r.snp
    #     shift_ind = ind
    #     if reference[shift_ind] != convert(DNA, r.ref)
    #         display((r, reference[shift_ind], convert(DNA, r.ref)))
    #         n += 1
    #     end
    # end
    # display(n)
    # return nrow(recurrent_df)
    # orfplot(recurrent_df)
    # top_n_snps = recurrent_df[1:10, :]
    # return format_ref_w_snps(top_n_snps.snp)
end
using BioSequences
function format_ref_w_snps(snp_list)
    reference_reader = FASTA.Reader(open(joinpath(@__DIR__, "../data/wuhCor1.fa")))
    reference = FASTA.sequence(only(reference_reader))[1:end-2]
    mutated_seq = deepcopy(reference)
    gene_starts = first.(sort(gene_index))
    for snp in snp_list
        (; ind, alt) = snp
        display(mutated_seq[ind])
        mutated_seq[ind] = alt
        display(mutated_seq[ind])
    end

    ref_aa = BioSequences.translate(reference)
    mutated_aa = BioSequences.translate(mutated_seq)# code=standard_genetic_code, allow_ambiguous_codons=true)
    for (; ind, alt) in snp_list
        gene_start = gene_starts[findlast(<=(ind), gene_starts)]
        aa_ind = floor(Int, ind / 3)
        display((reference[aa_ind*3:aa_ind*3+2], mutated_seq[aa_ind*3:aa_ind*3+2]))
        gene_aa_ind = floor(Int, (ind - gene_start) / 3)
        display((ref_aa[aa_ind], gene_aa_ind, mutated_aa[aa_ind]))
    end
    # do this biopython
    # grab coordinates
    # make vcf
    # seq_writer = FASTA.Writer(open(joinpath(@__DIR__, "../data/mutated_ref.fasta"), "w"))
    # rec = FASTA.Record("mutated", reference)
    # write(seq_writer, rec)
    # close(seq_writer)
    # close(reference_reader)
end


end
