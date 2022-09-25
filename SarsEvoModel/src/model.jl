
using LabelledArrays
using DifferentialEquations
using LoopVectorization
using Symbolics
using Polyester
using BenchmarkTools
using StaticArrays
import Base.length
@inline function sigma(x, y; sigma_x=4.0, sigma_y=4.0, rounding=true)
    if rounding
        round(exp(-1 * ((x)^2 / sigma_x + (y)^2 / sigma_y)); digits=1)
    else
        exp(-1 * ((x)^2 / sigma_x + (y)^2 / sigma_y))
    end
end


struct ModelParameters{T,T2,T3}
    begin_date::Date
    initial_population::Float64
    β::T
    ξ::Float64
    γ::Float64
    M::T2
    sigma_matrix::T3
    location_data::LocationData
    vaccination_matrix::Matrix{Float64}
    u0::Matrix{Float64}
end

function ModelParameters(
    begin_date,
    initial_population,
    β::T,
    ξ,
    γ,
    M::T2,
    sigma_matrix::T3,
    location_data,
    vaccination_matrix
) where {T,T2,T3}
    u0 = zeros(Float64, (w, h * 5))
    S = @view u0[1:w, 1:h]
    I = @view u0[1:w, (h+1):(h*2)]
    R = @view u0[1:w, (h*2+1):(h*3)]
    V = @view u0[1:w, (h*3+1):(h*4)]
    C = @view u0[1:w, (h*4+1):(h*5)]

    date_ind = findfirst(>=(begin_date), location_data.dates)
    stringency = location_data.stringency[date_ind:end]
    vaccination_mrna = location_data.vaccination_mrna[date_ind:end]
    cases_by_lineage = location_data.cases_by_lineage[date_ind:end]
    location_data_from_date = LocationData(
        location_data.dates[date_ind:end],
        stringency,
        cases_by_lineage,
        vaccination_mrna
    )
    init_mrna_vaccinated = sum(location_data.vaccination_mrna[1:date_ind])
    S .= initial_population
    mrna_vaccine = (x_i=2.7, y_i=3.8, width=20.0, pop=init_mrna_vaccinated)#B.1.1.7
    for (; x_i, y_i, width, pop) in (mrna_vaccine,), x in 1:w, y in 1:h
        x_i_transformed, y_i_transformed = map_coords_to_model_space(x_i, y_i)
        S[y, x] -= sigma(x - x_i_transformed, y - y_i_transformed; sigma_x=width, sigma_y=width, rounding=false) * pop
        V[y, x] += sigma(x - x_i_transformed, y - y_i_transformed; sigma_x=width, sigma_y=width, rounding=false) * pop
    end

    init_population_at_coords = [sum([m[x, y] for m in location_data.cases_by_lineage[1:date_ind]]) for x in 1:w, y in 1:h]
    active_population_at_coords = [sum([m[x, y] for m in location_data.cases_by_lineage[date_ind-6:date_ind]]) for x in 1:w, y in 1:h]
    I .+= active_population_at_coords
    S .-= (init_population_at_coords .+ active_population_at_coords)
    R += init_population_at_coords

    savefig(heatmap(I), plots_path("debug"))

    return ModelParameters{T,T2,T3}(
        begin_date,
        initial_population,
        β,
        ξ,
        γ,
        M,
        sigma_matrix,
        location_data_from_date,
        vaccination_matrix,
        u0,
    )
end
function Base.length(p::ModelParameters)
    return length(p.location_data.stringency)
end

function rhs(du, u, p, t, const_params)
    (; β, M, sigma_matrix) = p
    (; ξ, γ, initial_population, location_data, vaccination_matrix) = const_params
    (; stringency, vaccination_mrna) = location_data
    dS = @view du[1:w, 1:h]
    dI = @view du[1:w, (h+1):(h*2)]
    dR = @view du[1:w, (h*2+1):(h*3)]
    dV = @view du[1:w, (h*3+1):(h*4)]
    dC = @view du[1:w, (h*4+1):(h*5)]

    S = @view u[1:w, 1:h]
    I = @view u[1:w, (h+1):(h*2)]
    R = @view u[1:w, (h*2+1):(h*3)]
    hmone = h - 1
    wmone = w - 1
    day = trunc(Int, t)
    stringency_t = stringency[day+1]
    yesterday_vaccinations = day >= 1 ? vaccination_mrna[day] : 0.0
    vaccination_rate_by_day_t = 10 * (vaccination_mrna[day+1] - yesterday_vaccinations) / (initial_population * w * h)
    @inbounds for j in 2:hmone
        @inbounds for i in 2:wmone
            force_of_infection = zero(eltype(sigma_matrix))
            @turbo for l in 2:hmone
                for k in 2:wmone
                    force_of_infection += sigma_matrix[k-i+w, l-j+h] * I[k, l]
                end
            end
            dS[i, j] = -1 * stringency_t * β[i, j] * force_of_infection * S[i, j] + γ * R[i, j] - vaccination_rate_by_day_t * S[i, j]
            diffusion = M * (-4 * I[i, j] + I[i-1, j] + I[i+1, j] + I[i, j+1] + I[i, j-1])
            dI[i, j] = β[i, j] * stringency_t * I[i, j] * S[i, j] - ξ * I[i, j] + diffusion
            dR[i, j] = stringency_t * β[i, j] * (force_of_infection - I[i, j]) * S[i, j] + ξ * I[i, j] - γ * R[i, j]
            dV[i, j] = vaccination_rate_by_day_t * S[i, j] * vaccination_matrix[i, j]
            dC[i, j] = β[i, j] * stringency_t * I[i, j] * S[i, j] + diffusion
        end
    end
end

function import_callback(integrator, num_imports_per_day)
    omicron_x, omicron_y, width = antigenic_map_paper["B.1.1.529"]
    omicron_x, omicron_y = map_coords_to_model_space(omicron_x, omicron_y)
    S = @view integrator.u[1:w, 1:h]
    I = @view integrator.u[1:w, (h+1):(h*2)]
    for x in 1:w, y in 1:h
        S[y, x] -= sigma(x - omicron_x, y - omicron_y; sigma_x=width, sigma_y=width, rounding=false) * num_imports_per_day
        I[y, x] += sigma(x - omicron_x, y - omicron_y; sigma_x=width, sigma_y=width, rounding=false) * num_imports_per_day
    end
end



function create_model(params, const_params)
    #IC
    u0 = const_params.u0
    

    tstops = params.imports_start_time:(params.imports_start_time + 30.0)
    import_cb = DiscreteCallback((u, t, int) -> t in tstops, int -> import_callback(int, params.num_imports_per_day), save_positions=(false, false))
    f = ODEFunction((du, u, p, t) -> rhs(du, u, p, t, const_params))
    prob = ODEProblem(f, u0, (0.0, length(const_params) - 1), params)
    return prob, import_cb, tstops
end
