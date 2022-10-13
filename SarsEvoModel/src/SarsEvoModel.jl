module SarsEvoModel
export main
using BenchmarkTools
using Dates
using ProgressMeter
using Optimization, OptimizationBBO, StatsBase

using Setfield
using Arrow

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
include("data.jl")
include("model.jl")
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
const initial_pop = 329.5e6

function show_landscape(day)
    loc_data = USALocationData()
    p = [sum([m[x,y] for m in loc_data.cases_by_lineage]) for x in 1:size(loc_data.cases_by_lineage[1])[1], y in 1:size(loc_data.cases_by_lineage[1])[2]]
end


function main()

    begin_date = Date(2020, 9, 01)
    location_data = USALocationData()

    inv_infectious_period = 1 / 9
    r_0 = 2.1
    transmission_rate = r_0 * inv_infectious_period
    
    β = Float64[transmission_rate for i in 1:w, j in 1:h]
    x_mrna, y_mrna = map_coords_to_model_space(2.7, 3.8)
    const_params = ModelParameters(
        begin_date,
        initial_pop,
        β,
        inv_infectious_period,
        0.07,
        0.5,
        Float64[sigma(i, j) for i in -(w - 1):(w-1), j in -(h - 1):(h-1)],
        location_data,
        [sigma(i - x_mrna, j - y_mrna; sigma_x=20.0, sigma_y=20.0, rounding=false) for i in 1:w, j in 1:h]
    )
    incident_cases = sum.(const_params.location_data.cases_by_lineage)
    function optimization_objective(x)
        β = [x[1] + x[2] * i + x[3] * j for i in 1:w, j in 1:h] ./ initial_pop
        M = x[4]
        num_imports_per_day = x[7]
        imports_start_time = x[8]
        sigma_matrix = Float64[sigma(i, j; sigma_x=x[5], sigma_y=x[6]) for i in -(w - 1):(w-1), j in -(h - 1):(h-1)]
        prob, cb, tstops = create_model((; β, M, sigma_matrix, num_imports_per_day,imports_start_time), const_params)
        sol = DifferentialEquations.solve(prob, Tsit5(); callback=cb, saveat=1:1:length(const_params), tstops=tstops)
        yield()
        return sol
    end

    function loss(sol, x)
        if sol.retcode != :Success
            return Inf
        end
        err = 0
        umone = view(sol.u[1], 1:w, (h*4+1):(h*5))
        umone_sum = sum(umone)
        for (i, t) in enumerate(sol.t[2:end])
            ind = Int(t)
            u_i = view(sol.u[i+1], 1:w, (h*4+1):(h*5))
            u_i_sum = sum(umone_sum)

            # err += sum((const_params.location_data.cases_by_lineage[i] .- (u_i - umone)) .^ 2)
            err += (sum(incident_cases[ind]) - (u_i_sum - umone_sum))^2
            umone = u_i
        end
        err /= length(sol.t)
        yield()

        return err
    end

    x0 = [transmission_rate, 0.00, 0.00, 1.5, 4.0, 4.0, 20_000.0, 280.0]

    f = OptimizationFunction((x, _) -> loss(optimization_objective(x), x))
    prob = Optimization.OptimizationProblem(f, x0, 0; lb=[0.0, -0.5, -0.5, 0.01, 1.0, 1.0, 100.0,200.0], ub=[5.0, 0.5, 0.5, 5.0, 50.0, 50.0, 50_000.0,400.0])
    optimizers = ThreadsX.map(x -> Optimization.solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxtime=3000.0, verbose=true),1:Threads.nthreads())
    optimizer = argmin(o -> o.minimum, optimizers).u

    sol_opt = optimization_objective(optimizer)
    plot_solution(sol_opt, const_params)
    plot_parameters(optimizer)

    sol = optimization_objective(x0)
    plot_solution(sol, const_params)
    plot_parameters(x0)
    return optimizers
end

function plot_parameters(params)

    opt_β = [params[1] + params[2] * i + params[3] * j for i in 1:w, j in 1:h] ./ initial_pop
    plt = heatmap(opt_β; xlabel = "MDS1", ylabel = "MDS2", seriescolor = :Blues)
    savefig(plt, plots_path("model_beta"))

end
end