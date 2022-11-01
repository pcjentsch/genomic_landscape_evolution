module SarsEvoModel
export main
using BenchmarkTools
using Dates
using ProgressMeter
using Optimization, OptimizationBBO, StatsBase
using AverageShiftedHistograms

using Setfield
using Arrow


include("antigenic_mapping/antigenic_map.jl")
include("antigenic_mapping/genomic_types.jl")
include("data.jl")
include("plotting.jl")
include("antigenic_mapping/analyze_snps.jl")
include("antigenic_mapping/antigenic_map_plotting.jl")
include("antigenic_mapping/io.jl")

lower_triangular(n) = ((i, j) for i in 1:n, j in 1:n if j < i)
const datasets = (
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
        datapath("usa_sequences/alignments/all_lineages.csv"),
    )
)

function make_antigenic_map()
    binding_sites = Set(py"""list($(py_bcalc.sites))""") .+ S_gene_ind
    recurrent_df = parse_recurrent_mutations(datapath("filtered_mutations.tsv"))

    snp_weight_dict = make_snp_dict(recurrent_df, binding_sites)
    
    dataset = datasets[2]
        dataset_name, alignments, metadata, lineage_path = dataset

        @info "filtering unique genomes.."
        method = "homoplasy"
        genomes_w_metadata = serial_load_arrow(
            () -> get_data_fasta(alignments, metadata, binding_sites, lineage_path) ,
            datapath("cache","genomes_$dataset_name.arrow")
        )
        unique_df = serial_load_arrow(
            () -> unique_genomes(genomes_w_metadata, snp_weight_dict),
            datapath("cache","unique_df_$dataset_name.arrow")
        )
        filter_df =serial_load_arrow(
            () -> filter_genomes(unique_df),
            datapath("cache","filter_df_$dataset_name.arrow")
        )
        @info "computing pairwise distances $method"
        name = "$(method)_$(dataset_name)"
        dm = serial_load(
            () -> pairwise_distances(filter_df, snp_weight_dict),
            datapath("cache","$name.data")
        )
        manifold_projection(dm, filter_df)
        mds_df = serial_load_arrow(
                () -> innerjoin(
                select(filter_df,[:filtered_genome,:closest_mapped_lineage, :mds_x, :mds_y]),
                unique_df,
                on = [:filtered_genome, :closest_mapped_lineage]
            ),
            datapath("cache","unique_df_mds_$dataset_name.arrow")
        )
        # for k in 1:10
        #     filter_idxs = rand(1:nrow(filter_df), trunc(Int, nrow(filter_df) * 0.9))
        #     subsampled = filter_df[filter_idxs, :]
        #     subsampled_dm = dm[filter_idxs, filter_idxs]
        #     manifold_projection(subsampled_dm, subsampled)
        #     plot_mds("$name/subsample_$(k)_$name", subsampled, dm)
        # end
        # for lineage in unique(unique_df.closest_mapped_lineage)
        #     filter_idxs = filter_df.closest_mapped_lineage .!= lineage
        #     subsampled = filter_df[filter_idxs, :]
        #     subsampled_dm = dm[filter_idxs, filter_idxs]
        #     manifold_projection(subsampled_dm, subsampled)
        #     plot_mds("$name/filter_lineage_$(lineage)_$name", subsampled, dm)
        # end


        plot_mds("$name/$(name)_mds", mds_df)


        unique_df.week_submitted = round.(unique_df.Collection_Date, Week)
        bounds_x = extrema(unique_df.mds_x)
        bounds_y = extrema(unique_df.mds_y)
        x_grid = LinRange(bounds_x..., w)
        y_grid = LinRange(bounds_y..., h)
        anim = Animation()
        heatmap_anim = Animation()
        sort!(unique_df, :Collection_Date)
        @info "Interpolating..."
        grped_by_date = groupby(unique_df, :Collection_Date; sort=true)
        dates = Date[]
        grids = Vector{Matrix{Float64}}(undef, length(grped_by_date))


        @showprogress for (i, (key, gdf)) in enumerate(pairs(grped_by_date))
            o = ash(gdf.mds_x, gdf.mds_y; rngx=x_grid, rngy=y_grid, mx=10, my=10)
            grids[i] = o.z
            push!(dates, key.Collection_Date)
            p = plot()
            scatter!(p, gdf.mds_x, gdf.mds_y; xlims=bounds_x, ylims=bounds_y, markersize=1.5,
                markerstrokewidth=0.3,
                size=(400, 300),
                plotting_settings...)
            plot!(p, o.rngx, o.rngy, o.z; title=key.Collection_Date, xlims=bounds_x, ylims=bounds_y, plotting_settings...)
            htmp = heatmap(o.z; title=key.Collection_Date, plotting_settings...)
            frame(anim, p)
            frame(heatmap_anim, htmp)
        end
        gif(anim, plots_path("$name/$(name)_mds_density"; filetype="gif"))
        gif(heatmap_anim, plots_path("$name/$(name)_mds_density_heatmap"; filetype="gif"))
        serialize(datapath("cache","$(name)_grids.data"), (grids, dates))
end

function get_map_distance(lineage_a, lineage_b)
    mds_df = Arrow.Table(datapath("cache","unique_df_mds_usa.arrow"))
    ind_a = findfirst(==(lineage_a), mds_df.lineage)
    ind_b = findfirst(==(lineage_b), mds_df.lineage)
    (isnothing(ind_a) || isnothing(ind_b)) || error("lineages not present")
    return sqrt((mds_df.mds_x[ind_a] -mds_df.mds_x[ind_b])^2 + (mds_df.mds_y[ind_a] -mds_df.mds_y[ind_b])^2)
end
function stress_plot()
    homoplasy_dm = deserialize(datapath("cache","homoplasy_usa.data"))
    mstress = 8
    stress_list = zeros(mstress)
    p = plot()
    for (k, dm) in enumerate((homoplasy_dm,))
        display(k)
        for i in 1:mstress
            mds = fit(MDS, dm; distances=true, maxoutdim=i)
            stress_list[i] = stress(mds)
        end
        plot!(p, 1:mstress, stress_list; xlabel="MDS Out-dimension", ylabel="Stress", plotting_settings...)
    end
    savefig(p, plots_path("combined_mds_stress"))

end

end