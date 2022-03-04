module SarsEvoModel

using LabelledArrays
using DifferentialEquations
using LoopVectorization

@inline function sigma(x, y; sigma_x = 0.5, sigma_y = 0.5)
    return round(exp(-1 * ((x)^2 / sigma_x + (y)^2 / sigma_y)); digits = 1)
end

function rhs(du, u, p, t)
    (; β, ξ, γ, N, M, strain_dims, sigma_matrix) = p
    (xsize, ysize) = strain_dims .- 1

    dS = @view du[1:50, 1:50]
    dI = @view du[1:50, 51:100]
    dR = @view du[1:50, 101:150]
    S = @view u[1:50, 1:50]
    I = @view u[1:50, 51:100]
    R = @view u[1:50, 101:150]


    @inbounds for j in 2:ysize
        @inbounds for i in 2:xsize
            force_of_infection = zero(Float64)
            @inbounds for l in 2:ysize
                @inbounds for k in 2:xsize
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

using Symbolics
using BenchmarkTools
function main()
    strain_dims = (50, 50)



    u0 = zeros(Float64, strain_dims .* (1, 3))

    S = @view u0[1:50, 1:50]
    I = @view u0[1:50, 51:100]
    R = @view u0[1:50, 101:150]


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
    #unifrac
    β = (Float64).(LinRange(3.0, 8.0, 50))
    sigma_matrix = Float64[sigma(i - k, j - l) for i in axes(I, 1), j in axes(I, 2), k in axes(I, 1), l in axes(I, 2)]
    p = (β = β, ξ = 1.0, γ = 0.1, N = sum(u0) / length(S), sigma_matrix = sigma_matrix, M = 0.00001, strain_dims = strain_dims)

    du0 = copy(u0)
    rhs(du0, u0, p, 0.0)

    jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> rhs(du, u, p, 0.0), du0, u0)
    f = ODEFunction(rhs; jac_prototype = float.(jac_sparsity))
    prob = ODEProblem(f, u0, (0.0, 30_000.0), p)
    sol = solve(prob, Rodas5(); atol = 1e-10, rtol = 1e-10)

    return sol
end

using Plots
using BenchmarkTools

function test()

    p = heatmap(-15.0:15.0, -15.0:15.0, [sigma(i, j) for i in -15.0:15.0, j in -15.0:15.0])
    savefig(p, joinpath(@__DIR__, "../plots/$sigma.png"))
    sol = main()
    display(sol.t)

    anim = Animation()
    tlist = Float64[]
    S_ts = Float64[]
    I_ts = Float64[]
    R_ts = Float64[]
    max_t = maximum(sol.t)
    for (i, (u_t, t)) in enumerate(zip(sol.u, sol.t))
        S = @view u_t[1:50, 1:50]
        I = @view u_t[1:50, 51:100]
        R = @view u_t[1:50, 101:150]
        push!(S_ts, sum(S))
        push!(I_ts, sum(I))
        push!(R_ts, sum(R))
        push!(tlist, t)

        heatmap1 = heatmap(S, title = "S", seriescolor = cgrad(:Blues), clims = (0.0, 1.0))
        heatmap2 = heatmap(I, title = "I", seriescolor = cgrad(:Blues), clims = (0.0, 1.0))
        heatmap3 = heatmap(R, title = "R", seriescolor = cgrad(:Blues), clims = (0.0, 1.0))

        ts1 = plot(tlist, S_ts, title = "S", seriescolor = :Blue, xlims = (0.0, max_t))
        ts2 = plot(tlist, I_ts, title = "I", seriescolor = :Blue, xlims = (0.0, max_t))
        ts3 = plot(tlist, R_ts, title = "R", seriescolor = :Blue, xlims = (0.0, max_t))

        p = plot(heatmap1, heatmap2, heatmap3, ts1, ts2, ts3; layout = (2, 3), size = (1500, 600), dpi = 300)
        frame(anim, p)
    end
    gif(anim, joinpath(@__DIR__, "../plots/model.gif"))
end

end