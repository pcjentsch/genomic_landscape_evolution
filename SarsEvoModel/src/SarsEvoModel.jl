module SarsEvoModel
export main
using BenchmarkTools
using ProgressMeter
#Ontario population
const initial_pop = 329.5e6
#Size of antigenic grid
const w = 25
const h = 25
function map_coords_to_model_space(coords_x, coords_y; limx=(1, 10), limy=(1, 10))
    lx_l, lx_u = limx
    ly_l, ly_u = limy
    return trunc.(Int,
        (
            (coords_x - lx_l) / (lx_u - lx_l) * (w - 1) + 1,
            (coords_y - ly_l) / (ly_u - ly_l) * (h - 1) + 1,
        )
    )
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
        datapath("usa_sequences/usa_lineages.csv"),
    )
)

using Optimization
function main()

    third_wave_begin = Date(2021, 03, 01)
    third_wave_end = Date(2021, 08, 01)
    location_data = USALocationData()

    inv_infectious_period = 1 / 7
    r_0 = 2.1
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
    prob = create_model(location_data, third_wave_begin, third_wave_end, params)
    date_ind = findfirst(>=(third_wave_begin), location_data.dates)
    incident_cases = location_data.cases_by_lineage[date_ind:end]

    function optimization_objective(β)
        params = ModelParameters(
            β,
            1 / 10,
            0.07,
            initial_pop,
            0.00001,
            sigma_matrix,
        )
        new_p = (
            diffusion_kernel=prob.p.diffusion_kernel,
            params=params,
            stringency=prob.p.stringency,
            vaccination_mrna=prob.p.vaccination_mrna
        )

        new_prob = remake(prob; p=new_p)
        sol = solve(new_prob, Rodas5();)
        I_total = diff([sum(u_t[1:w, (h*4+1):(h*5)]) for u_t in sol.u])
        err = 0
        for (true_incident, incident) in zip(incident_cases, I_total)
            err += (true_incident - incident)^2
        end
        return err
    end
    # @show optimization_objective(rand(w, h))
    # @show optimization_objective(rand(w, h))

    sol = solve(prob, Rodas5(); saveat=1:length(incident_cases))
    plot_solution(sol, location_data, third_wave_begin)
    # plot_data(location_data)
    # # plot_antigenic_map()
    # return location_data
end



using AverageShiftedHistograms
function make_antigenic_map()
    binding_sites = Set(py"""list($(py_bcalc.sites))""") .+ S_gene_ind
    for dataset in datasets[2:2]
        name, alignments, metadata, lineage_path = dataset
        recurrent_df = parse_recurrent_mutations(datapath("filtered_mutations.tsv"))
        genomes_w_metadata = serial_load(
            () -> get_data_fasta(alignments, metadata, binding_sites, lineage_path),
            datapath("genomes_$name.data")
        )
        @info "filtering unique genomes.."
        unique_df, filter_df, snp_weight_dict = serial_load(
            () -> unique_genomes(genomes_w_metadata, recurrent_df, binding_sites),
            datapath("unique_df_$name.data")
        )
        @info "computing pairwise distances"

        unique_df, filter_df, dm = pairwise_distances(unique_df, filter_df, snp_weight_dict)
        @info "computing stress"

        plot_mds(name, unique_df, dm)

        unique_df.week_submitted = round.(unique_df.date_submitted, Week)

        bounds_x = extrema(unique_df.mds_x)
        bounds_y = extrema(unique_df.mds_y)
        x_grid = LinRange(bounds_x..., 25)
        y_grid = LinRange(bounds_y..., 25)
        anim = Animation()
        heatmap_anim = Animation()
        sort!(unique_df, :week_submitted)
        @info "Plotting..."
        grped_by_date = groupby(unique_df, :date_submitted; sort=true)
        dates = Date[]
        grids = Vector{Matrix{Float64}}(undef, length(grped_by_date))
        @showprogress for (i, (key, gdf)) in enumerate(pairs(grped_by_date))
            o = ash(gdf.mds_x, gdf.mds_y; rngx=x_grid, rngy=y_grid, mx=10, my=10)
            grids[i] = o.z
            push!(dates, key.date_submitted)
            p = plot()
            scatter!(p, gdf.mds_x, gdf.mds_y; xlims=bounds_x, ylims=bounds_y, markersize=1.5,
                markerstrokewidth=0.3,
                size=(400, 300),
                plotting_settings...)
            plot!(p, o.rngx, o.rngy, o.z; title=key.date_submitted, xlims=bounds_x, ylims=bounds_y, plotting_settings...)
            htmp = heatmap(o.z; title=key.date_submitted, plotting_settings...)
            frame(anim, p)
            frame(heatmap_anim, htmp)
        end
        gif(anim, plots_path("$(name)_mds_density"; filetype="gif"))
        gif(heatmap_anim, plots_path("$(name)_mds_density_heatmap"; filetype="gif"))
        serialize(SarsEvoModel.datapath("$(name)_grids.data"), (grids, dates))
        return unique_df, filter_df, grids
    end
end

# function map_snp_distances(filter_df)
#     limx = extrema(filter_df.mds_x)
#     limy = extrema(filter_df.mds_y)
#     filter_df.gridded_position = map((x, y) -> map_coords_to_model_space(x, y; limx, limy), filter_df.mds_x, filter_df.mds_y)
#     df = groupby(filter_df, :gridded_position)
#     return df
#     snp_distances_matrix = zeros(25, 25, 25, 25)
#     for ind in CartesianIndices(snp_distances_matrix)

#     end
# end




end
