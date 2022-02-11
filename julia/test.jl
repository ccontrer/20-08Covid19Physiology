
include("RiskCov19.jl")

using DataFrames
using OrdinaryDiffEq
using Plots, StatsPlots

import .RiskCov19

# DATA
data = RiskCov19.load_AHS_data([:Age], [:Num_comorbidities])
data = RiskCov19.RiskCov19Data(data.df[1:100, :], 
                               data.timecol, 
                               data.deadcol, 
                               data.lincovars, 
                               data.odecovars)
train, test = RiskCov19.data_partition(data, 0.7, shuffle=true)
# display
display(first(train.df, 6))
display(describe(train.df))
println("Dimensions: ", size(train.df))

# ODE Model
model = RiskCov19.model001
# display
sol = solve(model.prob, Tsit5(), dtmax=1e-1)
plot(sol, labels=model.vars.description[1], lw=2, ylims=(0., 1))

# Turing model
chn = RiskCov19.train(train, model, n_samples=300, n_chains=2)
# display
display(chn)
plot(chn)
# save results
savefig(plot(chn), joinpath("results", model.id*"-figure-chains.png"))
write(joinpath("results", model.id*"-chains.jls"), chn)
# predict
pred = RiskCov19.predict(train, chn, train, model)

# Turing model
chn = RiskCov19.train(train, model, n_samples=300, n_chains=2)
# display
display(chn)
plot(chn)
# save results
savefig(plot(chn), joinpath("results", model.id*"-figure-chains.png"))
write(joinpath("results", model.id*"-chains.jls"), chn)
# predict
pred = RiskCov19.predict(train, chn, train, model)