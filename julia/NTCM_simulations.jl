### Normal Tissue Comlication Model Simulations
### v0.1
##   - ODE model abstraction with its own structure (NTCM_ODE_models.jl)
##   - Logging into files (info and results)
### v0.2
##   - ODE solve our of patient-based loop
##
##   TODO: use MCMC to estimate parameters in the ODE (ForwardDiff issue, adjoint?)
### Notation
##   functions: lowerCammelNotation
##   variables: snake_notation


using OrdinaryDiffEq
using Turing
using LinearAlgebra
using QuadGK
using Plots, StatsPlots, LaTeXStrings
using DelimitedFiles, DataFrames

function hms(t)
  (h,r) = divrem(t, 60*60)
  (m,r) = divrem(r, 60)
  (s,r) = divrem(r, 1)
  string(Int(h), "h ", Int(m), "m ", Int(s), "s")
end


include("NTCM_ODE_models.jl")

figures_dir = "figures/";
results_dir = "results/"

## Model
model_ODE = model_VLF001
@info "MCMC simulation for model: $(model_ODE.id)
  log file: $(results_dir*model_ODE.id*"-logs.log")"
io = open(results_dir*model_ODE.id*"-logs.log", "w");

## Data
# load data
data, header = readdlm("NTCM_data.csv", ',', header=true)
df = DataFrame(data, vec(header))
df[!,:group] = convert.(Int64,df[:,:group])
df_groups_info = readdlm("NTCM_data_groups_info.csv", ',')[:]
write(io, "Data information: 
  number of entries: $(size(df, 1))
  number of columns: $(size(df, 2))
  number of groups: $(length(unique(df.group)))\n");
for (i, g) in enumerate(df_groups_info) 
  write(io, "\t$(i) => $(g)\n");
end

## Model
write(io, "ODE model info:
  Model id: $(model_ODE.id)
  Model name: $(model_ODE.name)
  Variables: $(model_ODE.vars[:labels])
  Covariate variables: $(model_ODE.vars[:labels][model_ODE.vars[:covars]])
  Number of ODE parameters: $(length(model_ODE.prob.p))
    Fixed parameters ($(length(model_ODE.pars[:fixed][:labels]))): $(model_ODE.pars[:fixed][:labels])
    Parameters with a prior ($(length(model_ODE.pars[:prior][:labels]))): $(model_ODE.pars[:prior][:labels])
    Randomized parameters ($(length(model_ODE.pars[:random][:labels]))): $(model_ODE.pars[:random][:labels])\n");

# Probabilistic model
write(io, "Probabilistic model:
  Covariates: E (healthy lung cells)
  Fitting parameters: β (covariates coeficients)
  Fitting parameters (groups): λ (hazard rate)
  Randomized parameters: t_s (latent infectious time)");
close(io);

@model function modelProb(data, model_ODE::NTCM_model)
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
    pf = pars[:fixed][:values]
    pp = pars[:prior][:values]
    # doesn't work
    # η = Vector{Float64}(undef, length(pp))
    # for i = 1:length(pp)
    #   η[i] ~ truncated(Normal(pp[i], 0.001), 0.0, Inf) # doesn't work
    # end
    # η ~ filldist(Exponential(3.7e3), length(pp))
    η = pp
    pr = pars[:random][:values]
    ζ = Vector(undef, length(pr))
    for i = 1:length(pr)
      ζ[i] = rand(truncated(Normal(pr[i], 0.001), 0.0, Inf))
    end
    # ν random effects (groups)
    # ν ~ 
    ## Solve for each group not each patient
    # latent time
    pθ = pars[:thresholds][:values]
    t_s = rand(Normal(5.0, 0.1))
    # ODE
    prob = remake(model_ODE.prob, p=vcat(pf, η, ζ), tspan=(0.0, t_s + maximum(data.time)))
    sol = solve(prob, Tsit5(), dtmax=1e-1)
    
    # Covariate parameters
    covars = model_ODE.vars[:covars]
    N_covars = length(covars)
    μ_β = fill(0., N_covars)
    Σ_β = 0.1*I
    β ~ MvNormal(μ_β, Σ_β)
    
    # iterate for each patient
    for row in eachrow(data)
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
n_samples = 1000
n_chains = 4
sampler = NUTS()
model = modelProb(df, model_ODE)
tstat = @timed chns = mapreduce(c -> sample(model, sampler, n_samples), chainscat, 1:n_chains)

@info "Simulation completed (runtime: $(hms(tstat.time)))
  chain save in file: $(results_dir*model_ODE.id*"-chains.log")"
savefig(plot(chns), results_dir*model_ODE.id*"-figure-chains.png")
write(results_dir*model_ODE.id*"-chains.jls", chns)
