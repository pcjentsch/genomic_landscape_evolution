module SarsEvoModel
export main
using BenchmarkTools

#Ontario population
const initial_pop = 14.57e6
#Size of antigenic grid
const w = 50
const h = 50
function map_coords_to_model_space(coords_x, coords_y)
    return trunc.(Int, (coords_x * (w / 10), coords_y * (h / 10)))
end
include("antigenic_map.jl")
include("model.jl")
include("plotting.jl")
include("data.jl")




function main()
    init_data = SarsEvoModel.load_covid_data(Date(2021, 03, 01), Date(2021, 07, 01))
    # β = Float64[0.015 * round(exp(-1 * ((i - 25)^2 / 1e3 + (j - 25)^2 / 1e3)), digits=2) for i in 1:w, j in 1:h]
    # sigma_matrix = Float64[sigma(i - k, j - l) for i in 1:w, j in 1:h, k in 1:w, l in 1:h]

    # params = ModelParameters(
    #     β,
    #     0.01,
    #     0.07,
    #     initial_pop,
    #     0.00001,
    #     sigma_matrix,
    # )
    # sol = run(init_data, Date(2021, 03, 01), params)
    # plot_solution(sol)

    # plot_antigenic_map()
end


end