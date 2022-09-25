# # Logistic regression

# First, we import DynamicHMC and related libraries,

using TransformVariables, LogDensityProblems, DynamicHMC, TransformedLogDensities, MultivariateStats

# then some packages that help code the log posterior,

using Parameters, Statistics, Random, Distributions, LinearAlgebra, StatsFuns, LogExpFunctions

# then diagnostic and benchmark tools,

using MCMCDiagnosticTools, BenchmarkTools, ManifoldLearning

# and use ForwardDiff for AD since the dimensions is small.

import ForwardDiff



struct bMDS{Tp,Tδ}
    n::Int
    d::Matrix{Tδ}
    p::Tp
    function bMDS(d::Matrix{Tδ}, p::Tp) where {Tδ, Tp}
        return new{Tp,Tδ}(size(d)[1],d,p)
    end
end
function normcdf(x)
    cdf(Normal(), x)
end
function (problem::bMDS)(θ)
    (;n, d, p) = problem
    (;xs,σ) = θ
    ll = 0.0
    # for i in 1:n
    #     ll += loglikelihood(MvNormal(zeros(p), I(p)), xs[i])
    # end

    ll += loglikelihood(Exponential(),p)

    m = div(n*(n-1),2)
    ssr = 0.0
    lognorm_term = 0.0
    for i in 1:n
        for j in 1:(i-1)
            δ_ij = sqrt(sum((xs[i] .- xs[j]).^2))
            ssr += (d[i,j] - δ_ij)^2
            lognorm_term += log(normcdf(δ_ij/σ))
        end
    end
    return (σ^2)^(-m/2) * exp(-1*(1/(2*σ^2)) * ssr - lognorm_term)
end

#test dataset
test_dataset = ManifoldLearning.scurve()[1]
dist_matrix(pts) = reshape([sum((pts[:,i] .- pts[:,j]).^2) for i in 1:size(pts)[2], j in 1:size(pts)[2]],(size(pts)[2],size(pts)[2]))
dm = dist_matrix(test_dataset)
bmds = bMDS(dm, 2)


# mds = fit(MDS,dm; maxoutdim = 2, distances = true)
# pts = predict(mds)
 
# scatter(pts[1,:], pts[2,:])

for i in 1:10
    test_xs = [rand(mds.p) for i in 1:mds.n]
    test_sigma = 1.0
    test_p = 1.0
    @show mds((xs = test_xs, σ = test_sigma, p = test_p))
end
# # Make up parameters, generate data using random draws.

# N = 1000
# X = hcat(ones(N), randn(N))
# y = rand.(Bernoulli.(logistic.(X*β)));

# # Create a problem, apply a transformation, then use automatic differentiation.

# p = MDS(dm, 2)   # data and (vague) priors
# t = as((xs = as(Array, length(test_xs)), σ = as(Real, -∞,∞), p = as(Real, -∞,∞) )) # identity transformation, just to get the dimension
# P = TransformedLogDensity(t, p)      # transformed
# ∇P = ADgradient(:ForwardDiff, P)

# # Sample using NUTS, random starting point.

# results = map(_ -> mcmc_with_warmup(Random.default_rng(), ∇P, 1000), 1:5)

# # Extract the posterior. (Here the transformation was not really necessary).

# β_posterior = first.(transform.(t, eachcol(pool_posterior_matrices(results))))

# # Check that we recover the parameters.

# mean(β_posterior)

# # Quantiles

# qs = [0.05, 0.25, 0.5, 0.75, 0.95]
# quantile(first.(β_posterior), qs)

# quantile(last.(β_posterior), qs)

# # Check that mixing is good.

# ess, R̂ = ess_rhat(stack_posterior_matrices(results))