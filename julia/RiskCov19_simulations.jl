### Risk of Covid-19 infections - Simulations
### v0.1
##   - ODE model abstraction with its own structure (NTCM_ODE_models.jl)
##   - Logging into files (info and results)
### v0.2
##   - ODE solve our of patient-based loop
### v0.3
##   - 
##
##   TODO: use MCMC to estimate parameters in the ODE (ForwardDiff issue, adjoint?)
### Notation
##   functions: UowerCammelNotation
##   variables: snake_notation


using Distributions
using OrdinaryDiffEq
using Turing
using LinearAlgebra
using QuadGK
using Plots, StatsPlots, LaTeXStrings 
using DelimitedFiles, DataFrames
using CategoricalArrays

function hms(t)
  (h,r) = divrem(t, 60*60)
  (m,r) = divrem(r, 60)
  (s,r) = divrem(r, 1)
  string(Int(h), "h ", Int(m), "m ", Int(s), "s")
end

figures_dir = "figures/"
results_dir = "results/"

## Model
println("Loading models")
include("RiskCov19_ODE_models.jl")
model_ODE = model_VLF001
println("Using model '$(model_ODE.name)'")
@info "MCMC simulation for model: $(model_ODE.id)
  log file: $(results_dir*model_ODE.id*"-logs.log")"
io = open(results_dir*model_ODE.id*"-logs.log", "w");

## Survival model
include("RiskCov19_survival_models.jl")

## Data
# load data
data_filename = pwd()*"/../../data/AHS/Restricted/analysis.csv"
data, header = readdlm(data_filename, ',', header=true)
df = DataFrame(data, vec(header))
# handeling data types
df = identity.(df) # auto identify the data type
features = ["Disch_days", "Dead", 
            "Age", "Age group", "Sex", 
            "Chronic Pulmonary Disease", "Num. comorbidities"]
select!(df, features)
# Remove missing
deleteat!(df, df.Sex.=="")
# Transform
rename!(df, replace.(names(df), " " => "_"))
rename!(df, replace.(names(df), "." => ""))
transform!(df, :Disch_days => ByRow(x -> Float64(-x)) => :Time,
               :Dead => ByRow(x -> x == "True" ? 1.0 : 0.0) => :Dead,
               :Sex => ByRow(x -> x == "Male" ? 1.0 : 0.0) => :Male,
               :Chronic_Pulmonary_Disease => ByRow(x -> x == "True" ? 1 : 0),
               :Num_comorbidities => ByRow(x -> Float64(x));
               renamecols=false)
select!(df, Not([:Disch_days, :Sex]))
df."Age_group" = categorical(df."Age_group", ordered=true)
levels!(df."Age_group", ["Under 1 year", "1-4 years", "5-9 years", "10-19 years", "20-29 years", 
        "30-39 years", "40-49 years", "50-59 years", "60-69 years", "70-79 years", "80+ years"])
onehot!(df, :Age_group)
# 
data = RiskCov19Data(df, :Time, :Dead, [:Male, :Sex, :No_comorbidities], [:Chronic_Pulmonary_Disease])
# data info
nrows, ncols = size(data.df)
write(io, "Data information: 
  number of entries: $(nrows)
  number of features: $(ncols - 2)\n");
write(io, "  Covariates\n");
for str in data.covarcols
  write(io, "\t$(str)\n");
end
write(io, "  Ode groups\n");
for str in data.odegrcols
  write(io, "\t$(str)\n");
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


@model function modelProb(data::RiskCov19Data, ODEModel::RiskCov19ODEModel, ::Type{T} = Float64) where {T}
    # data: (time, dead_flag, groups)
    # group: K columns
    # covars: indexes in model_ODE.prob.u of covariates eg. [1, ]
    
    # groups
    df = data.df
    N = size(df, 1) 
    t = df[:, data.timecol]
    d = df[:, data.deadcol]
    X = Matrix(df[:, data.covarcols])
    num_covars = length(data.covarcols)
    Y = Matrix(df[:, data.odegrcols])
    odegrs = unique(Y)
    num_odegrs = length(odegrs)
    
    # Basal hazard rate
    λ ~ Exponential(1.0)
    if λ <= 0
        Turing.@addlogprob! -Inf
        return
    end
    hazard_dist = Exponential(λ)

    # Covariates
    # β ~ filldist(Normal(0., 1), num_covars)
    μ_β = fill(0., num_covars)
    Σ_β = 1 ./ std(X, dims=1)[:]
    β ~ MvNormal(μ_β, Σ_β)
    
    # # ODE parameters (default + random sample)
    # pars = ODEModel.pars
    # pf = pars[:fixed][:values]
    # pp = pars[:prior][:values]
    # pr = pars[:random][:values]
    # pθ = pars[:thresholds][:values]
    # η = Vector{Vector{T}}(undef, num_odegrs)
    # sol = Vector{ODESolution}(undef, num_odegrs)
    # # latent time
    t_s = rand(Normal(5.0, 0.1))
    # for k in 1:num_odegrs
    #   η[k] = Vector{T}(undef, length(pp))
    #   for i = 1:length(pp)
    #     η[k][i] ~ truncated(Normal(pp[i], 0.001), 0.0, Inf)
    #   end
    #   # ODE
    #   prob = remake(ODEModel.prob, p=vcat(pf, η[k], pr), tspan=(0.0, t_s + maximum(t)))
    #   sol[k] = solve(prob, Tsit5())
    # end
    # # Y = reshape(1.0 .- map(t -> sol(t)[1], t), :, 1)
    # γ ~ filldist(Normal(0., 0.1), 1)
    
    x = X * β
    
    # iterate for each patient
    for i in 1:N
        # user Log-likelihood
        # z = (1. - sol[Y[i]+1](t[i])[1]) * γ
        Turing.@addlogprob! logL(t_s + t[i], d[i], x[i], hazard_dist)
    end
end

## Simulation
n_samples = 1000
n_chains = 4
sampler = NUTS()
model = modelProb(data, model_ODE)
tstat = @timed chns = mapreduce(c -> sample(model, sampler, n_samples), chainscat, 1:n_chains)

plot(chns[:,:,1])

@info "Simulation completed (runtime: $(hms(tstat.time)))
  chain save in file: $(results_dir*model_ODE.id*"-chains.log")"
savefig(plot(chns), results_dir*model_ODE.id*"-figure-chains.png")
write(results_dir*model_ODE.id*"-chains.jls", chns)
