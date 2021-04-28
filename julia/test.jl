using Revise

using Optim, LsqFit
using Plots, StatsPlots, LaTeXStrings
using CSV
using Statistics
using Turing
using DataFrames
using Distributions

df = DataFrame(CSV.File("./R/dat.csv"))

f(x, a, b, c) = a + (b - a)*exp(-exp(c)*x)

@model function model(data)
     ~ Normal(25, sqrt(2))
    ϕ₁ ~ Normal(8, sqrt(2))
    ϕ₂ ~ Normal(-1, sqrt(2))

    γ₀ ~ Normal(0, 1)
    γ₁ ~ Normal(1, 1)
    γ₂ ~ Normal(0, 1)
    
    τϕ₁ ~ TDist(3)
    τϕ₂ ~ TDist(3)
    τϕ₃ ~ TDist(3)
    τγ₁ ~ TDist(3)
    τγ₂ ~ TDist(3)
    τγ₃ ~ TDist(3)

    Ω ~ LKJ(1)

    for i ∈ eachindex(data)
        for t ∈  eachindex(datai)
            σ[i,t] = exp(f(trial, γ₀, γ₁, γ₂))
            ϵ[i,t] ~ Normal(0, σ[i,j])
            y[i,t] = f(trial, ϕ₀, ϕ₁, ϕ₂) + ϵ[i, t]
        end
    end
end


@model function test(x)
    a ~ Normal(0, 1)
    b ~ Normal(5, 1)
    # c ~ Normal(0, 1) + Uniform(0, 1)
    d = a + b
    x = d
end

chain = sample(test(missing), Prior(), 1000)