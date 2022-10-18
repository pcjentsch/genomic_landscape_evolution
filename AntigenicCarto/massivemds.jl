using RCall, ManifoldLearning, Plots, LinearAlgebra, Serialization, InlineStrings, SarsEvoModel
test_dataset = first(ManifoldLearning.scurve(1000))

dist_matrix(pts) = reshape([sum((pts[:,i] .- pts[:,j]).^2) for i in 1:size(pts)[2], j in 1:size(pts)[2]],(size(pts)[2],size(pts)[2]))

# dm = dist_matrix(test_dataset)
# locs = rand(1000,2)
burnin = 1_000
nsamples = 1_000

R"devtools::load_all(\"~/Work/MassiveMDS\")"
(_, filter_df, dm) = deserialize("../SarsEvoModel/data/homoplasy_usa.data")
small_dm = dm[1:1000,1:1000]
locs = permutedims(mapreduce(x -> [x...],hcat, filter_df.mapped_lineage_position[1:1000])) .+ (rand(1000,2) .- 0.5)
@rput small_dm
@rput locs
#this works!!
R"""
rlocs <- locs
hmc <- MassiveMDS::hmcsampler(n_iter=$nsamples +$burnin, data=small_dm, burnIn=$burnin, learnPrec=TRUE, learnTraitPrec=TRUE, gpu=1, treeCov=FALSE, trajectory=0.1, locations = rlocs)
"""
@rget hmc 
samples = Vector{Matrix{Float64}}(hmc[:samples])

# function logprior(X, U, V, Uinv, Vinv)
#     n,p = size(X)

#     grad_prior = -0.5 * (transpose(Vinv * transpose(X) * Uinv) + Uinv * X * Vinv)
#     product = Vinv * transpose(X) * Uinv * X
#     prior_exponent = -0.5 * sum(diag(product))

#     logdet_U = logdet(U)
#     logdet_V = logdet(V)

#     prior = prior_exponent - (n*p / 2) * log(2*pi) - (n/2) * logdet_V - (p/2) * logdet_U

#     return prior, grad_prior
# end


# #lets try it with HMC coded in julia now
# using AdvancedHMC, Distributions

# const U = Diagonal(ones(1000))
# const Uinv = Diagonal(ones(1000))
# const V = Diagonal(ones(2))
# const Vinv = Diagonal(ones(2))
# @rput U
# @rput Uinv
# @rput V
# @rput Vinv
# function log_likelihood_grad(θ)
#     reshaped_theta=reshape(θ, (1000,2))
#     R"""
#     rθ <- $reshaped_theta
#     engine <- MassiveMDS::updateLocations(engine, rθ)
#     log_likelihood <- MassiveMDS::getLogLikelihood(engine)
#     log_likelihoodgrad <- MassiveMDS::getGradient(engine)
#     """
#     @rget log_likelihood
#     @rget log_likelihoodgrad
#     priormat, gradpriormat = logprior(reshaped_theta, U,V, Uinv, Vinv)
#     return -1 .* (priormat + log_likelihood), vec(-(gradpriormat + log_likelihoodgrad))
# end
# function log_likelihood(θ)
#     reshaped_theta=reshape(θ, (1000,2))
#     R"""
#     rθ <- $reshaped_theta
#     engine <- MassiveMDS::updateLocations(engine, rθ)
#     log_likelihood <- MassiveMDS::getLogLikelihood(engine)
#     """
#     @rget log_likelihood
#     priormat, gradpriormat = logprior(reshaped_theta, U,V, Uinv, Vinv)
#     return -1 .* (priormat + log_likelihood)
# end



# R"""
# rlocs <- locs
# engine <- MassiveMDS::createEngine(2, 1000, TRUE, 8, 2, 1, FALSE)
# engine <- MassiveMDS::setPairwiseData(engine, as.matrix(dm))
# engine <- MassiveMDS::updateLocations(engine, rlocs)
# engine <- MassiveMDS::setPrecision(engine, precision = 1)
# """
# # Choose parameter dimensionality and initial parameter value
# D = (1000*2); initial_θ = rand(D)

# # Define the target distribution
# # ℓπ(θ) = logpdf(MvNormal(zeros(D), I), θ)

# # Set the number of samples to draw and warmup iterations
# n_samples, n_adapts = 5_000, 2_000

# # Define a Hamiltonian system
# metric = DiagEuclideanMetric(D)
# hamiltonian = Hamiltonian(metric, log_likelihood, log_likelihood_grad)

# # Define a leapfrog solver, with initial step size chosen heuristically
# initial_ϵ = find_good_stepsize(hamiltonian, initial_θ)
# integrator = Leapfrog(initial_ϵ)

# # Define an HMC sampler, with the following components
# #   - multinomial sampling scheme,
# #   - generalised No-U-Turn criteria, and
# #   - windowed adaption for step-size and diagonal mass matrix
# proposal = NUTS{MultinomialTS, GeneralisedNoUTurn}(integrator)
# adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(0.8, integrator))

# Run the sampler to draw samples from the specified Gaussian, where
#   - `samples` will store the samples
#   - `stats` will store diagnostic statistics for each sample
# samples, stats = sample(hamiltonian, proposal, initial_θ, n_samples, adaptor, n_adapts; progress=true)
# p = plot()
# for sample in samples
#     reshaped = reshape(sample, (1000,2))
#     scatter!(p, reshaped[:,1],reshaped[:,2],markersize = 0.1,markerstrokewidth=0.0,legend=false )
# end
# display(p)