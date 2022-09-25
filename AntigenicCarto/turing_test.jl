using Turing, StatsPlots, Random, LinearAlgebra, MultivariateStats, ManifoldLearning, ReverseDiff
using Plots
using BenchmarkTools
Turing.setadbackend(:reversediff)

@model function bmds(p,d, ::Type{T} = Float64) where {T}
    n = size(dm)[1]
    xs = Vector{Vector{T}}(undef,n)
    #priors for variance
    σ2 ~ Exponential()
    #priors for points, just put them at zero with wide prior
    for i in 1:n
        xs[i] ~ MvNormal(zeros(p),Diagonal(fill(10,p)))
    end

    for i in 1:n
        for j in 1:(i-1)
            δ_ij = 0.0
            for k in 1:p
                δ_ij += (xs[i][k] - xs[j][k])^2
            end
            d[i,j] ~ truncated(Normal(sqrt(δ_ij),σ2); lower = 0.0)
        end
    end
end

@model function bmds_faster(p,d, ::Type{T} = Float64) where {T}
    n = size(dm)[1]
    m = div(n*(n-1),2)
    xs = Vector{Vector{T}}(undef,n)
    #priors for variance
    σ2 ~ Exponential()
    #priors for points, just put them at zero with wide prior
    for i in 1:n
        xs[i] ~ MvNormal(zeros(p),Diagonal(fill(10,p)))
    end
    δ = zeros(T, m)
    ind = 1
    for i in 1:n
        for j in 1:(i-1)
            δ_ij = zero(T)
            for k in 1:p
                δ_ij  += (xs[i][k] - xs[j][k])^2
            end
            δ[ind] = exp(sqrt(δ_ij)) #exp transform
            ind += 1
        end
    end
    d ~ MvLogNormal(δ,Diagonal(fill(σ2,m)))
end

function make_test_dataset()
    pts = ManifoldLearning.scurve()[1]
    
    n = size(pts)[2]
    m = div(n*(n-1),2)
    dists = Vector{Float64}(undef,m)
    ind = 1
    for i in 1:n
        for j in 1:(i-1)
            dists[ind] = sqrt(sum((pts[:,i] .- pts[:,j]).^2))
            ind +=1
        end
    end
    return dists
end
#test dataset

dm = make_test_dataset()
chain = sample( bmds_faster(2, dm), HMC(0.1,5), MCMCThreads(), 10, 4)
# describe(chain)

# # Plot a summary of the sampling process for the parameter p, i.e. the probability of heads in a coin.
# histogram(chain[:p])
    