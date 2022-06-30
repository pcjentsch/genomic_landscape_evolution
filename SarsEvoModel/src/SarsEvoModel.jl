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
using Optimization, ForwardDiff, ReverseDiff, OptimizationNLopt, OptimizationBBO, StatsBase
function main()

    begin_date = Date(2020, 9, 01)
    location_data = USALocationData()

    inv_infectious_period = 1 / 7
    r_0 = 1.55
    transmission_rate = r_0 * inv_infectious_period
    diffusion_kernel = [sigma(i, j; sigma_x=20.0, sigma_y=20.0, rounding=false) for i = -25:25, j = -25:25]
    diffusion_kernel = diffusion_kernel ./ sum(diffusion_kernel)
    β = Float64[transmission_rate for i in 1:w, j in 1:h]
    sigma_matrix = Float64[sigma(i - k, j - l) for i in 1:w, j in 1:h, k in 1:w, l in 1:h]
    initial_pop = 329.5e6

    params = ModelParameters(
        begin_date,
        initial_pop,
        β,
        inv_infectious_period,
        0.07,
        0.0001,
        sigma_matrix,
        location_data,
        diffusion_kernel
    )
    u0 = params.u0
    du0 = similar(u0, size(u0))

    jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> rhs(du, u, params, 0.0), du0, u0)
    incident_cases = sum.(params.location_data.cases_by_lineage)
    # display(plot(incident_cases))
    function optimization_objective(x)
        new_β = reshape(x, (w, h))
        new_p = @set params.β = new_β
        prob = create_model(new_p, jac_sparsity)
        sol = solve(prob, Rodas5(); saveat=1:length(params))
        return sol
    end

    function loss(sol, x)
        if sol.retcode != :Success
            return Inf
        end
        β = reshape(x, (w, h))
        l = length(sol.u)
        err = 0
        # reconstructed = Float64[]
        umone = sum(sol.u[1][1:w, (h*4+1):(h*5)])
        for i in 2:l

            u_i = sum(sol.u[i][1:w, (h*4+1):(h*5)])
            # push!(reconstructed, (u_i - umone))

            err += (incident_cases[i] - (u_i - umone))^2
            umone = u_i
        end
        # p = plot(reconstructed; ylims=(0, max(maximum(reconstructed), maximum(incident_cases))))

        # plot!(incident_cases; ylims=(0, max(maximum(reconstructed), maximum(incident_cases))))
        # display(p)
        err /= l
        yield()
        # display(err)
        # display(total_variation(β) * 1e6)
        return err + total_variation(β) * 1e5
    end
    # @show ForwardDiff.derivative(x -> optimization_objective(x, 0), 0.0001)
    f = OptimizationFunction((x, _) -> loss(optimization_objective(x), x))
    β = Float64[transmission_rate for i in 1:w, j in 1:h]
    x0 = vec(β)
    @show loss(optimization_objective(x0), x0)

    prob = Optimization.OptimizationProblem(f, x0, 0, lb=fill(0.0001, length(x0)), ub=fill(transmission_rate * 1.2, length(x0)))
    optimizer = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters=100000, maxtime=5000.0)
    sol = optimization_objective(optimizer)

    # # sol = optimization_objective(x0)
    # # return loss(sol, x0)
    plot_solution(sol, params)
    # loss(sol, x0)
    # plot_data(location_data)
    # # plot_antigenic_map()
    # return location_data
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
