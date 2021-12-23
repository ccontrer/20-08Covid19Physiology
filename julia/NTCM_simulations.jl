### Normal Tissue Comlication Model Simulations
### v0.1
## Notation
##   functions: lowerCammelNotation
##   variables: snake_notation

using OrdinaryDiffEq
using Turing
using LinearAlgebra
using QuadGK
using Plots, StatsPlots, LaTeXStrings
using DelimitedFiles, DataFrames

include("NTCM_ODE_models.jl")

figures_dir = "figures/";
results_dir = "results/"

## Data
# load data
data, header = readdlm("NTCM_data.csv", ',', header=true)
df = DataFrame(data, vec(header))
df[!,:group] = convert.(Int64,df[:,:group])
# df = DataFrame(CSV.File("NTCM_data.csv")) # using CSV
df_groups_info = readdlm("NTCM_data_groups_info.csv", ',')[:]
@info "Data information: " rows=size(df, 1) cols=size(df, 2) groups=length(unique(df.group))

## ODE model
model_ODE = model_ODE_VLF001
@info "ODE model info:
  number of variables: $(length(model_ODE.prob.u0))
  number of parameters: $(length(model_ODE.prob.p))" vars=u_labels

# Probabilistic model
@info "Probabilistic model:
  Covariates: E (healthy lung cells)
  Fitting parameters: β (covariates coeficients)
  Fitting parameters (groups): λ (hazard rate)
  Randomized parameters: t_s (latent infectious time)
  Randomized parameters (groups): μ_E (birth rate of healthy lung cells)"
@model function modelProb(data, model_ODE, covars)
    # data: (time, dead_flag, group)
    # group: value 1 to K_groups 
    # covars: indexes in model_ODE.prob.u of covariates eg. [1, ]
    
    # groups
    K_groups = length(unique(data.group))
    
    # ODE parameters (default)
    
    # priors
    #Vmax ~ Normal(3.7e3, 10.)
    # μ_E = rand(Normal(model_ODE.pars[:prior][1], 0.0001), K_groups)
    λ ~ filldist(Exponential(1.0), K_groups)
    
    # ODE parameters (default + random sample)
    pars = model_ODE.pars
    pf = collect(values(pars[:fixed]))
    pp = collect(values(pars[:prior]))
    # doesn't work
    # η = Vector{Float64}(undef, length(pp))
    # for i = 1:length(pp)
    #   η[i] ~ truncated(Normal(pp[i], 0.001), 0.0, Inf) # doesn't work
    # end
    # η ~ filldist(Exponential(3.7e3), length(pp))
    η = pp
    pr = collect(values(pars[:random]))
    ζ = Vector(undef, length(pr))
    for i = 1:length(pr)
      ζ[i] = rand(truncated(Normal(pr[i], 0.001), 0.0, Inf))
    end
    # ν random effects (groups)
    # ν ~ 
    pθ = collect(values(pars[:thresholds]))
    
    # Covariate parameters
    N_covars = length(covars)
    μ_β = fill(0., N_covars)
    Σ_β = 0.1*I
    β ~ MvNormal(μ_β, Σ_β)
    
    # iterate for each patient
    for row in eachrow(data)
        # latent time
        t_s = rand(Normal(5.0, 0.1))
        # ODE
        prob = remake(model_ODE.prob, p=vcat(pf, η, ζ), tspan=(0.0, t_s + row.time))
        sol = solve(prob, Tsit5(), dtmax=1e-1)
        # Survival
        x(t) = sol(t)[covars, :] - pθ
        h(t) = λ[row.group]*exp(-dot(x(t), β))
        H(t) = quadgk(h, 0., t)
        # Log-likelihood
        logL(t, d) = sum(d.*log(h(t)) .- H(t))
        Turing.@addlogprob! logL(t_s + row.time, row.dead)
    end
end


## Simulation
n_samples = 50
n_chains = 1
model = modelProb(df, model_ODE, [1, ])
chns = sample(model, NUTS(), n_samples)
# chns = sample(model, NUTS(), MCMCThreads(), n_samples, n_chains)
plot(chns)
savefig(plot(chns), results_dir*"figure-chains-"*model_ODE.id)
write(results_dir*"chains-"*model_ODE.id*".jls", chns)
