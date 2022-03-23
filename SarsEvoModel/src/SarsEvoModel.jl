module SarsEvoModel
export main
using BenchmarkTools
include("antigenic_map.jl")
include("model.jl")
include("plotting.jl")
include("data.jl")

function main()
    init_data = SarsEvoModel.load_covid_data(Date(2021, 03, 01), Date(2021, 07, 01))

    sol = run(init_data, Date(2021, 03, 01))
    plot_solution(sol)

    plot_antigenic_map()
end


end