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
include("model.jl")
include("data.jl")
include("plotting.jl")
include("antigenic_mapping/analyze_snps.jl")
include("antigenic_mapping/antigenic_map_plotting.jl")
include("antigenic_mapping/genomic_types.jl")
include("antigenic_mapping/io.jl")



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


function make_antigenic_map()

    datasets = (
        ("uk", joinpath(@__DIR__, "../data/uk_sequences/alignments/"), joinpath(@__DIR__, "../data/uk_sequences/uk_sequences_metadata_new.tsv")),
        ("usa", joinpath(@__DIR__, "..data/usa_sequences/alignments/"), joinpath(@__DIR__, "../data/usa_sequences/usa_sequences_metadata.tsv"))
    )
    binding_sites = Set(py"""list($(py_bcalc.sites))""") .+ S_gene_ind
    # name, alignments, metadata = datasets[1]
    # genomes_w_metadata = get_data_fasta(alignments, metadata, binding_sites)
    recurrent_df = parse_recurrent_mutations(joinpath(@__DIR__, "../data/filtered_mutations.tsv"))
    orfplot(recurrent_df)
    return recurrent_df[1:10, :]
    # plot_samples(genomes_w_metadata)
    # genomes_w_metadata_new, unique_genomes_df_new, mds = analyze_binding_retained(genomes_w_metadata, recurrent_df)
    # lineages_df = load_lineages()
    # plot_mds(genomes_w_metadata_new, unique_genomes_df_new, lineages_df)
    # end
end


end
