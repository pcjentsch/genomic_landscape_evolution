
using LabelledArrays
using DifferentialEquations
using LoopVectorization
using Symbolics
using BenchmarkTools
using StaticArrays

@inline function sigma(x, y; sigma_x = 2.5, sigma_y = 2.5)
    return round(exp(-1 * ((x)^2 / sigma_x + (y)^2 / sigma_y)); digits = 1)
end

function rhs(du, u, p, t)
    (; β, ξ, γ, N, M, sigma_matrix) = p

    dS = @view du[1:w, 1:h]
    dI = @view du[1:w, (h+1):(h*2)]
    dR = @view du[1:w, (h*2+1):(h*3)]
    S = @view u[1:w, 1:h]
    I = @view u[1:w, (h+1):(h*2)]
    R = @view u[1:w, (h*2+1):(h*3)]


    @inbounds for j in 2:(h-1)
        @inbounds for i in 2:(w-1)
            force_of_infection = zero(Float64)
            @inbounds for l in 2:(h-1)
                @inbounds for k in 2:(w-1)
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
const w = 50
const h = 50

function run(init_data)

    u0 = zeros(Float64, (w, h * 3))

    S = @view u0[1:w, 1:h]
    I = @view u0[1:w, (h+1):(h*2)]
    R = @view u0[1:w, (h*2+1):(h*3)]


    #IC
    I[25, 25] = 0.2
    S .= 0.9
    S[1, :] .= 0
    S[end, :] .= 0
    S[:, 1] .= 0
    S[:, end] .= 0

    I[1, :] .= 0.0
    I[end, :] .= 0.0
    I[:, 1] .= 0.0
    I[:, end] .= 0.0

    for (lineage, (x0, y0, width), population) in init_data
        for x in 1:w, y in 1:h
            I[x, y] += sigma(x - x0, y - y0; sigma_x = width, sigma_y = width) .* population
        end
    end

    return I
    #unifrac
    # β = SVector{50}((Float64).(LinRange(4.0, 4.0, 50)))
    # sigma_matrix = Float64[sigma(i - k, j - l) for i in axes(I, 1), j in axes(I, 2), k in axes(I, 1), l in axes(I, 2)]
    # p = (β = β, ξ = 1.0, γ = 0.1, N = sum(u0) / length(S), sigma_matrix = sigma_matrix, M = 0.001, strain_dims = strain_dims)

    # du0 = copy(u0)
    # # @btime rhs($du0, $u0, $p, 0.0)

    # jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> rhs(du, u, p, 0.0), du0, u0)

    # f = ODEFunction(rhs; jac_prototype = float.(jac_sparsity))

    # prob = ODEProblem(f, u0, (0.0, 200.0), p)

    # sol = solve(prob, Rosenbrock23();
    #     # atol = 1e-10,
    #     # rtol = 1e-10,
    #     progress = true
    # )

    # return sol
end
