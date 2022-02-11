struct RiskCov19ODEVars
    names::Vector{Symbol}
    description::Vector{String}
    ICs::Vector{AbstractFloat}
    covars_ix::Vector{Integer}
end

struct RiskCov19ODEPars
    names::Vector{Symbol}
    values::Vector{AbstractFloat}
end

struct RiskCov19ODEModel
    id::String
    name::String
    prob::ODEProblem
    vars::RiskCov19ODEVars
    pars::RiskCov19ODEPars
end

## Virus load function
function VLF(t::Float64, p::Vector)
    a₁, a₂, b₁, b₂, α, Vmax, Vmin = p
    H = 0.5*(tanh(1e3*(t - a₂)) + 1.)
    v₁ = 1.0 + (Vmax - 1.0)*(tanh(6.0*(t - (a₁ + a₂)/2)/(a₂ - a₁)) - tanh(-3.0*(a₂ + a₁)/(a₂ - a₁)))/2
    v₂ = 1.0 - H + H*exp(-α*(t - a₂))
    v₃ = 1.0 - (1.0 - Vmin)*(tanh(6.0*(t - (b₁ + b₂)/2)/(b₂ - b₁)) - tanh(-3.0*(b₂ + b₁)/(b₂ - b₁)))/2
    v₁*v₂*v₃
end

## Model: VLF001
## VLF + Lung tissue
# ODE model
function modelODE!(du, u, p, t)
    E, = u
    a₁ = 0.5     # a₁
    a₂ = 4.0     # a₂
    b₁ = 13.0    # b₁
    b₂ = 19.0    # b₂
    α = 0.1     # α
    Vmin = 1e-8    # Vmin
    Vmax = 3.7e3   # Vmax
    μ_E, δ_E = p
    V(t) = VLF(t, [a₁, a₂, b₁, b₂, α, Vmax, Vmin])
    # Lung damage
    du[1] = dE = μ_E*E*(1.0 - E) - δ_E*E*V(t)/Vmax # Healthy lung epithelium
end
u0 = [
    1.0   #E: Healthy lung epithilium
]
vars = RiskCov19ODEVars([:E, ], ["Healthy lung cells"], u0, [1, ])
tspan = (0.0, 40.0)
pars = RiskCov19ODEPars([:μ_E, :δ_E], [0.01, 0.15])
prob = ODEProblem(modelODE!, u0, tspan, pars.values)
model001 = RiskCov19ODEModel("VLF001", "VLF + lung tissue", prob, vars, pars)