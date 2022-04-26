module SarsEvoModel
export main
using BenchmarkTools

#Ontario population
const initial_pop = 67.22e6
#Size of antigenic grid
const w = 25
const h = 25
function map_coords_to_model_space(coords_x, coords_y)
    return trunc.(Int, (coords_x * (w / 10), coords_y * (h / 10)))
end
include("antigenic_map.jl")
include("model.jl")
include("data.jl")
include("plotting.jl")




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


end