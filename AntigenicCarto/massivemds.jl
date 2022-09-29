using RCall, ManifoldLearning
test_dataset = ManifoldLearning.scurve(1000)[1]
dist_matrix(pts) = reshape([sum((pts[:,i] .- pts[:,j]).^2) for i in 1:size(pts)[2], j in 1:size(pts)[2]],(size(pts)[2],size(pts)[2]))
dm = dist_matrix(test_dataset)
locs = rand(1000,2)
@rput dm
@rput locs
#this works!!
# R"""
# rlocs <- locs
# hmc <- MassiveMDS::hmcsampler(n_iter=1000, data=dm, burnIn=0, learnPrec=TRUE, learnTraitPrec=TRUE, gpu=1, treeCov=FALSE, trajectory=0.1, locations = rlocs)
# """
# @rget hmc 
# samples = Vector{Matrix{Float64}}(hmc[:samples])


#lets try it with HMC coded in julia now
using AdvancedHMC, Distributions
using LinearAlgebra


R"""
rlocs <- locs
engine <- MassiveMDS::createEngine(2, 1000, TRUE, 8, 2, 0, FALSE)
engine <- MassiveMDS::setPairwiseData(engine, as.matrix(dm))
engine <- MassiveMDS::updateLocations(engine, rlocs)
engine <- MassiveMDS::setPrecision(engine, precision = 1)
"""
function log_likelihood_grad(θ)
    reshaped_theta=reshape(θ, (1000,2))

    R"""
    rθ <- $reshaped_theta
    engine <- MassiveMDS::updateLocations(engine, rθ)
    log_likelihood <- MassiveMDS::getLogLikelihood(engine)
    log_likelihoodgrad <- MassiveMDS::getGradient(engine)
    """
    @rget log_likelihood
    @rget log_likelihoodgrad
    return log_likelihood, vec(log_likelihoodgrad)
end
function log_likelihood(θ)
    reshaped_theta=reshape(θ, (1000,2))
    R"""
    rθ <- $reshaped_theta
    engine <- MassiveMDS::updateLocations(engine, rθ)
    log_likelihood <- MassiveMDS::getLogLikelihood(engine)
    """
    @rget log_likelihood
    return log_likelihood
end


# Choose parameter dimensionality and initial parameter value
D = (1000*2); initial_θ = rand(D)

# Define the target distribution
# ℓπ(θ) = logpdf(MvNormal(zeros(D), I), θ)

# Set the number of samples to draw and warmup iterations
n_samples, n_adapts = 2_000, 1_000

# Define a Hamiltonian system
metric = DiagEuclideanMetric(D)
hamiltonian = Hamiltonian(metric, log_likelihood, log_likelihood_grad)

# Define a leapfrog solver, with initial step size chosen heuristically
initial_ϵ = find_good_stepsize(hamiltonian, initial_θ)
integrator = Leapfrog(initial_ϵ)

# Define an HMC sampler, with the following components
#   - multinomial sampling scheme,
#   - generalised No-U-Turn criteria, and
#   - windowed adaption for step-size and diagonal mass matrix
proposal = NUTS{MultinomialTS, GeneralisedNoUTurn}(integrator)
adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(0.8, integrator))

# Run the sampler to draw samples from the specified Gaussian, where
#   - `samples` will store the samples
#   - `stats` will store diagnostic statistics for each sample
samples, stats = sample(hamiltonian, proposal, initial_θ, n_samples, adaptor, n_adapts; progress=true)