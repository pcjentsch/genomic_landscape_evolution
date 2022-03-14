
using LabelledArrays
using DifferentialEquations
using LoopVectorization
using Symbolics
using BenchmarkTools
using StaticArrays

const w = 50
const h = 50
@inline function sigma(x, y; sigma_x = 2.5, sigma_y = 2.5, rounding = true)
    if rounding
        round(exp(-1 * ((x)^2 / sigma_x + (y)^2 / sigma_y)); digits = 1)
    else
        exp(-1 * ((x)^2 / sigma_x + (y)^2 / sigma_y))
    end
end

function rhs(du, u, p, t)
    (; β, ξ, γ, N, M, sigma_matrix) = p

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
                    force_of_infection += sigma_matrix[k, l, i, j] * I[k, l]
                end
            end
            dS[i, j] = -(β[j] / N) * force_of_infection * S[i, j] + γ * R[i, j]
            dI[i, j] = (β[j] / N) * I[i, j] * S[i, j] - ξ * I[i, j] +
                       M * (-4 * I[i, j] + I[i-1, j] + I[i+1, j] + I[i, j+1] + I[i, j-1])
            dR[i, j] = ξ * I[i, j] - γ * R[i, j]
        end
    end
end

function run(init_data, begin_date)

    u0 = zeros(Float64, (w, h * 3))

    S = @view u0[1:w, 1:h]
    I = @view u0[1:w, (h+1):(h*2)]
    R = @view u0[1:w, (h*2+1):(h*3)]

    #IC
    S .= 14.57e6

    for (lineage, (x0, y0, width), population_by_date) in init_data
        init_population_at_coords = sum(filter(:date => <=(begin_date), population_by_date).pop)
        active_population_at_coords = sum(filter(:date => x -> begin_date <= x <= begin_date + Day(5), population_by_date).pop)
        for x in 1:w, y in 1:h
            I[y, x] += sigma(x - x0, y - y0; sigma_x = width, sigma_y = width, rounding = false) .* active_population_at_coords
            S[y, x] -= sigma(x - x0, y - y0; sigma_x = width, sigma_y = width, rounding = false) .* init_population_at_coords + sigma(x - x0, y - y0; sigma_x = width, sigma_y = width) .* active_population_at_coords
            R[y, x] += sigma(x - x0, y - y0; sigma_x = width, sigma_y = width, rounding = false) .* init_population_at_coords
        end
    end

    # return S, I, R
    β = (Float64).(LinRange(0.1, 0.1, h))
    sigma_matrix = Float64[sigma(i - k, j - l) for i in 1:w, j in 1:h, k in 1:w, l in 1:h]
    p = (β = β, ξ = 0.01,
        γ = 0.05, N = sum(S) / length(S),
        sigma_matrix = sigma_matrix, M = 0.0001)

    du0 = zeros(size(u0))
    @btime rhs($du0, $u0, $p, 0.0)
    jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> rhs(du, u, p, 0.0), du0, u0)
    f = ODEFunction(rhs; jac_prototype = float.(jac_sparsity))

    prob = ODEProblem(f, u0, (0.0, 500.0), p)

    sol = solve(prob, Rosenbrock23();
        progress = true
    )

    return sol
end
