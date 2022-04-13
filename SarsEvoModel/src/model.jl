
using LabelledArrays
using DifferentialEquations
using LoopVectorization
using Symbolics
using BenchmarkTools
using StaticArrays

@inline function sigma(x, y; sigma_x=50.0, sigma_y=50.0, rounding=true)
    if rounding
        round(exp(-1 * ((x)^2 / sigma_x + (y)^2 / sigma_y)); digits=1)
    else
        exp(-1 * ((x)^2 / sigma_x + (y)^2 / sigma_y))
    end
end

struct ModelParameters
    β::Matrix{Float64}
    ξ::Float64
    γ::Float64
    initial_population::Float64
    M::Float64
    sigma_matrix::Array{Float64,4}
end


function rhs(du, u, p, t)
    (; β, ξ, γ, initial_population, M, sigma_matrix) = p

    dS = @view du[1:w, 1:h]
    dI = @view du[1:w, (h+1):(h*2)]
    dR = @view du[1:w, (h*2+1):(h*3)]
    S = @view u[1:w, 1:h]
    I = @view u[1:w, (h+1):(h*2)]
    R = @view u[1:w, (h*2+1):(h*3)]
    hmone = h - 1
    wmone = w - 1

    @inbounds for j in 2:hmone
        @inbounds for i in 2:wmone
            force_of_infection = zero(eltype(sigma_matrix))
            @inbounds for l in 2:hmone
                @inbounds for k in 2:wmone
                    force_of_infection += (β[i, j] / initial_population) * sigma_matrix[k, l, i, j] * I[k, l]
                end
            end
            dS[i, j] = -1 * force_of_infection * S[i, j] + γ * R[i, j]
            dI[i, j] = (β[i, j] / initial_population) * I[i, j] * S[i, j] - ξ * I[i, j] +
                       M * (-4 * I[i, j] + I[i-1, j] + I[i+1, j] + I[i, j+1] + I[i, j-1])
            dR[i, j] = ξ * I[i, j] - γ * R[i, j]
        end
    end
end

function run(location_data, begin_date, end_date, params::ModelParameters)
    (; β, sigma_matrix, initial_population) = params

    u0 = zeros(Float64, (w, h * 3))

    S = @view u0[1:w, 1:h]
    I = @view u0[1:w, (h+1):(h*2)]
    R = @view u0[1:w, (h*2+1):(h*3)]

    #IC
    S .= initial_population
    date_ind = findfirst(>=(begin_date), location_data.dates)

    init_population_at_coords = sum(location_data.cases_by_lineage[1:date_ind])
    active_population_at_coords = sum(location_data.cases_by_lineage[1:date_ind+5])
    I .+= active_population_at_coords
    S .-= (init_population_at_coords .+ active_population_at_coords)
    R += init_population_at_coords

    display(init_population_at_coords)
    # TODO make this better
    #identify strain with closest neutralization and pick that as centre for gaussian for vaccination IC
    # init_mrna_vaccinated = initial_population * 0.3 #sum(filter(:date => <=(begin_date), population_by_date).pop)
    # init_az_vaccinated = initial_population * 0.1 #sum(filter(:date => <=(begin_date), population_by_date).pop)
    # mrna_vaccine = (x_i=2.7, y_i=3.8, width=30.0, pop=init_mrna_vaccinated)#B.1.1.7
    # az_vaccine = (x_i=2.7, y_i=3.8, width=20.0, pop=init_az_vaccinated)
    # for (; x_i, y_i, width, pop) in (mrna_vaccine, az_vaccine), x in 1:w, y in 1:h
    #     x_i_transformed, y_i_transformed = map_coords_to_model_space(x_i, y_i)
    #     S[y, x] -= sigma(x - x_i_transformed, y - y_i_transformed * 5; sigma_x=width, sigma_y=width, rounding=false) * pop
    #     R[y, x] += sigma(x - x_i_transformed, y - y_i_transformed * 5; sigma_x=width, sigma_y=width, rounding=false) * pop
    # end


    ##Plot sigma, beta, and S_0 for debugging
    plt1 = heatmap(-(w // 2):w//2, -h//2:h//2, sigma_matrix[:, :, 25, 25], xlabel="antigenic distance", ylabel="antigenic distance"; plot_options...)
    savefig(plt1, plots_path("sigma"))
    plt2 = heatmap(1:w, 1:h, β; plot_options...)
    savefig(plt2, plots_path("beta"))
    plt3 = heatmap(1:w, 1:h, S; plot_options...)
    savefig(plt3, plots_path("S_0"))

    # du0 = zeros(size(u0))
    # jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> rhs(du, u, params, 0.0), du0, u0)
    # f = ODEFunction(rhs; jac_prototype=float.(jac_sparsity))

    # prob = ODEProblem(f, u0, (0.0, 60_000.0), params)

    # sol = solve(prob, Rodas5();
    #     progress=true
    # )
    # return sol
end
