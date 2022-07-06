module SarsEvoModel
export main
using BenchmarkTools
using Dates
using ProgressMeter
using Setfield
#Ontario population
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
        datapath("usa_sequences/usa_lineages.csv"),
    )
)
function total_variation(A)
    t = 0.0
    for i in 1:first(size(A))-1, j in 1:last(size(A))-1
        t += abs(A[i+1, j] - A[i, j]) + abs(A[i, j] - A[i, j+1])
    end
    return t
end
using Optimization, ForwardDiff, ReverseDiff, OptimizationNLopt, OptimizationBBO, StatsBase, Zygote
using SciMLSensitivity, ProfileView


function main()

    begin_date = Date(2020, 9, 01)
    location_data = USALocationData()

    inv_infectious_period = 1 / 7
    r_0 = 1.65
    transmission_rate = r_0 * inv_infectious_period
    initial_pop = 329.5e6
    β = Float64[transmission_rate for i in 1:w, j in 1:h]
    x_mrna, y_mrna = map_coords_to_model_space(2.7, 3.8)
    const_params = ModelParameters(
        begin_date,
        initial_pop,
        β,
        inv_infectious_period,
        0.07,
        0.5,
        Float64[sigma(x, y; sigma_x=4, sigma_y=4) for x in -(w - 1):(w-1), y in -(h - 1):(h-1)],
        location_data,
        [sigma(i - x_mrna, j - y_mrna; sigma_x=20.0, sigma_y=20.0, rounding=false) for i in 1:w, j in 1:h]
    )
    u0 = const_params.u0
    du0 = similar(u0, size(u0))

    # jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> rhs(du, u, params, 0.0), du0, u0)
    # jac_sparsity = float.(Symbolics.jacobian_sparsity((du, u) -> rhs(du, u, β, 0.0, const_params), du0, u0))
    incident_cases = sum.(const_params.location_data.cases_by_lineage)
    # display(plot(incident_cases))
    function optimization_objective(x)
        β = [x[1] + x[2] * i + x[3] * j for i in 1:w, j in 1:h] ./ const_params.initial_population
        M = x[4]
        sigma_matrix = Float64[sigma(i, j; sigma_x=x[5], sigma_y=x[6]) for i in -(w - 1):(w-1), j in -(h - 1):(h-1)]
        prob = create_model((; β, M, sigma_matrix), const_params)
        sol = solve(prob, Tsit5(); saveat=1:1:length(const_params)) #every 7 days
        return sol
    end

    function loss(sol, x)
        if sol.retcode != :Success
            return Inf
        end
        err = 0
        umone = view(sol.u[1], 1:w, (h*4+1):(h*5))
        umone_sum = sum(umone)
        # a = Animation()
        for (i, t) in enumerate(sol.t[2:end])
            ind = Int(t)
            u_i = view(sol.u[i+1], 1:w, (h*4+1):(h*5))
            u_i_sum = sum(umone_sum)
            # p = heatmap(const_params.location_data.cases_by_lineage[i])
            # p2 = heatmap(u_i .- umone)
            # fig = plot(p, p2)
            # display(fig)
            # frame(a, fig)
            err += sum((const_params.location_data.cases_by_lineage[i] .- (u_i - umone)) .^ 2)
            # err += (sum(incident_cases[ind]) - (u_i_sum - umone_sum))
            umone = u_i
        end
        # gif(a, plots_path("debug"; filetype=".gif"))
        err /= length(sol.t)
        yield()
        return err
    end
    x0 = [transmission_rate, -0.001, 0.000, 0.01, 4.0, 4.0]
    # @btime $loss($optimization_objective($x0), $x0)
    # β = [x0[1] + x0[2] * i + x0[3] * j for i in 1:w, j in 1:h] ./ const_params.initial_population
    # M = x0[4]

    # @btime rhs($du0, $u0, $((; β, M, const_params.sigma_matrix)), 0.0, $const_params)
    # @profview foreach(x -> rhs(du0, u0, ((; β, M, const_params.sigma_matrix)), 0.0, const_params), 1:10_000)

    f = OptimizationFunction((x, _) -> loss(optimization_objective(x), x))
    prob = Optimization.OptimizationProblem(f, x0, 0; lb=[0.001, -0.5, -0.5, 0.01, 0.1, 0.1], ub=[1.0, 0.5, 0.5, 5.0, 40.0, 40.0])#lb=fill(0.001, length(x0)), ub=fill(transmission_rate * 1.2, length(x0)))

    optimizers = ThreadsX.map(
        x -> solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters=1000000, maxtime=100000.0)
        Threads.nthreads + 1)

    return optimizers
    # sol = optimization_objective(optimizer)
    # plot_solution(sol, const_params)
    # loss(sol, const_params)
end

#2522.98 secs, 2777 evals, 2657 steps, improv/step: 0.250 (last = 0.0000), fitness=222894924.461761177
# ^C┌ Warning: Optimization interrupted, recovering intermediate results...
# └ @ BlackBoxOptim /home/workaccount/.julia/packages/BlackBoxOptim/dvDGl/src/opt_controller.jl:483
# 588
# Progress: 100%|██████████████████████████████████████████████| Time: 0:01:42
# ┌ Info: Saved animation to
# └   fn = "/home/workaccount/Work/postdoc/VoC_model/SarsEvoModel/plots/model.gif"
# u: 6-element Vector{Float64}:
#   0.48816792944503556
#  -0.03378768040205102
#  -0.007711016471248729
#   0.3594928335707491
#   0.5209967948289511
#  18.46451481226305

end