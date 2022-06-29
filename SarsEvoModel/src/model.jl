
using LabelledArrays
using DifferentialEquations
using LoopVectorization
using Symbolics
using BenchmarkTools
using StaticArrays

@inline function sigma(x, y; sigma_x=4.0, sigma_y=4.0, rounding=true)
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
    (; diffusion_kernel, params, stringency, vaccination_mrna) = p
    (; β, ξ, γ, initial_population, M, sigma_matrix) = params
    dS = @view du[1:w, 1:h]
    dI = @view du[1:w, (h+1):(h*2)]
    dR = @view du[1:w, (h*2+1):(h*3)]
    dV = @view du[1:w, (h*3+1):(h*4)]
    dC = @view du[1:w, (h*4+1):(h*5)]

    S = @view u[1:w, 1:h]
    I = @view u[1:w, (h+1):(h*2)]
    R = @view u[1:w, (h*2+1):(h*3)]
    V = @view u[1:w, (h*3+1):(h*4)]
    hmone = h - 1
    wmone = w - 1
    day = trunc(Int, t)
    stringency_t = stringency[day+1]
    yesterday_vaccinations = day >= 1 ? vaccination_mrna[day] : 0.0
    vaccination_rate_by_day_t = (vaccination_mrna[day+1] - yesterday_vaccinations) / initial_population
    @inbounds for j in 2:hmone
        @inbounds for i in 2:wmone
            force_of_infection = zero(eltype(sigma_matrix))
            @inbounds for l in 2:hmone
                @inbounds for k in 2:wmone
                    force_of_infection += (β[i, j] / initial_population) * sigma_matrix[k, l, i, j] * I[k, l]
                end
            end

            dS[i, j] = -1 * stringency_t * force_of_infection * S[i, j] + γ * R[i, j] - vaccination_rate_by_day_t * S[i, j]
            diffusion = M * (-4 * I[i, j] + I[i-1, j] + I[i+1, j] + I[i, j+1] + I[i, j-1])
            dI[i, j] = (β[i, j] / initial_population) * stringency_t * I[i, j] * S[i, j] - ξ * I[i, j] + diffusion
            dR[i, j] = ξ * I[i, j] - γ * R[i, j]
            x_i_transformed, y_i_transformed = map_coords_to_model_space(2.7, 3.8)
            dV[i, j] = vaccination_rate_by_day_t * S[i, j] * sigma(i - x_i_transformed, j - y_i_transformed; sigma_x=20.0, sigma_y=20.0, rounding=false)
            dC[i, j] = (β[i, j] / initial_population) * stringency_t * I[i, j] * S[i, j] + diffusion

        end
    end
end

# function rhs_diffusion_kernel(du, u, p, t)
#     (; diffusion_kernel, params, stringency, vaccination_mrna) = p
#     (; β, ξ, γ, initial_population, M, sigma_matrix) = params
#     dS = @view du[1:w, 1:h]
#     dI = @view du[1:w, (h+1):(h*2)]
#     dR = @view du[1:w, (h*2+1):(h*3)]
#     dV = @view du[1:w, (h*3+1):(h*4)]
#     dC = @view du[1:w, (h*4+1):(h*5)]

#     S = @view u[1:w, 1:h]
#     I = @view u[1:w, (h+1):(h*2)]
#     R = @view u[1:w, (h*2+1):(h*3)]

#     hmone = h - 1
#     wmone = w - 1
#     day = trunc(Int, t)
#     stringency_t = stringency[day+1]
#     yesterday_vaccinations = day >= 1 ? vaccination_mrna[day] : 0.0
#     vaccination_rate_by_day_t = 0.1 * (vaccination_mrna[day+1] - yesterday_vaccinations) / (initial_population * w * h)
#     for j in 2:hmone
#         for i in 2:wmone
#             force_of_infection = zero(eltype(sigma_matrix))
#             for l in 2:hmone
#                 for k in 2:wmone
#                     force_of_infection += (β[i, j] / initial_population) * sigma_matrix[k, l, i, j] * I[k, l]
#                 end
#             end
#             diffusion = 0.0
#             for l in 2:hmone
#                 for k in 2:wmone
#                     diffusion += I[k, l] * diffusion_kernel[i-k+wmone, j-l+hmone] # sigma(i - k, j - l; sigma_x=20.0, sigma_y=20.0, rounding=false)
#                 end
#             end

#             dS[i, j] = -1 * stringency_t * force_of_infection * S[i, j] + γ * R[i, j] - vaccination_rate_by_day_t * S[i, j]
#             dI[i, j] = (β[i, j] / initial_population) * stringency_t * I[i, j] * S[i, j] - ξ * I[i, j] + M * diffusion
#             dR[i, j] = ξ * I[i, j] - γ * R[i, j]

#             x_i_transformed, y_i_transformed = map_coords_to_model_space(2.7, 3.8)
#             dV[i, j] = vaccination_rate_by_day_t * S[i, j] * sigma(i - x_i_transformed, j - y_i_transformed; sigma_x=20.0, sigma_y=20.0, rounding=false)
#             dC[i, j] = M * diffusion + (β[i, j] / initial_population) * stringency_t * I[i, j] * S[i, j]
#         end
#     end

# end
using OffsetArrays
function create_model(location_data, begin_date, end_date, params)
    (; β, sigma_matrix, initial_population) = params
    (; dates, stringency, vaccination_mrna) = location_data
    u0 = zeros(Float64, (w, h * 5))

    S = @view u0[1:w, 1:h]
    I = @view u0[1:w, (h+1):(h*2)]
    R = @view u0[1:w, (h*2+1):(h*3)]
    V = @view u0[1:w, (h*3+1):(h*4)]
    C = @view u0[1:w, (h*4+1):(h*5)]
    #IC
    S .= initial_population
    date_ind = findfirst(>=(begin_date), location_data.dates)

    init_population_at_coords = [sum([m[x, y] for m in location_data.cases_by_lineage[1:date_ind]]) for x in 1:w, y in 1:h]
    active_population_at_coords = [sum([m[x, y] for m in location_data.cases_by_lineage[date_ind:date_ind+5]]) for x in 1:w, y in 1:h]
    I .+= active_population_at_coords
    display(sum(I))

    S .-= (init_population_at_coords .+ active_population_at_coords)
    R += init_population_at_coords
    stringency = stringency[date_ind:end]
    vaccination_mrna = vaccination_mrna[date_ind:end]
    init_mrna_vaccinated = sum(vaccination_mrna[1:date_ind])
    mrna_vaccine = (x_i=2.7, y_i=3.8, width=20.0, pop=init_mrna_vaccinated)#B.1.1.7
    # az_vaccine = (x_i=2.7, y_i=3.8, width=20.0, pop=init_az_vaccinated)
    for (; x_i, y_i, width, pop) in (mrna_vaccine,), x in 1:w, y in 1:h
        x_i_transformed, y_i_transformed = map_coords_to_model_space(x_i, y_i)
        S[y, x] -= sigma(x - x_i_transformed, y - y_i_transformed * 5; sigma_x=width, sigma_y=width, rounding=false) * pop
        V[y, x] += sigma(x - x_i_transformed, y - y_i_transformed * 5; sigma_x=width, sigma_y=width, rounding=false) * pop
    end

    diffusion_kernel = [sigma(i, j; sigma_x=20.0, sigma_y=20.0, rounding=false) for i = -25:25, j = -25:25]
    diffusion_kernel = diffusion_kernel ./ sum(diffusion_kernel)

    ##Plot sigma, beta, and S_0 for debugging
    # plt1 = heatmap(-div(w, 2):div(w, 2), -div(h, 2):div(h, 2), sigma_matrix[:, :, div(w, 2), div(h, 2)], xlabel="antigenic distance", ylabel="antigenic distance"; plotting_settings...)
    # savefig(plt1, plots_path("sigma"))
    # plt2 = heatmap(1:w, 1:h, β; plotting_settings...)
    # savefig(plt2, plots_path("beta"))
    # plt3 = heatmap(1:w, 1:h, S; plotting_settings...)
    # savefig(plt3, plots_path("S_0"))

    du0 = zeros(size(u0))
    jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> rhs(du, u, (; diffusion_kernel, params, stringency, vaccination_mrna), 0.0), du0, u0)
    f = ODEFunction(rhs; jac_prototype=float.(jac_sparsity))
    prob = ODEProblem(f, u0, (0.0, length(stringency) - 1), (; diffusion_kernel, params, stringency, vaccination_mrna))
    return prob
end
