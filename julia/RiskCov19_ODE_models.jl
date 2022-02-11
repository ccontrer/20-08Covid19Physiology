### Normal Tissue Comlication Model ODE models
### v0.1

struct RiskCov19ODEModel
    id::String
    name::String
    prob::ODEProblem
    vars::Dict
    pars::Dict
end

function parsDict2Array(pars::Dict)
  vcat(pars[:fixed][:values], pars[:prior][:values], pars[:random][:values])
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
    # p = [pf, pp, pr]
    E, = u
    a₁, a₂, b₁, b₂, α, Vmin, Vmax, μ_E, δ_E = p
    V(t) = VLF(t, [a₁, a₂, b₁, b₂, α, Vmax, Vmin])
    # Lung damage
    du[1] = dE = μ_E*E*(1.0 - E) - δ_E*E*V(t)/Vmax # Healthy lung epithelium
end
u0 = [
    1.0   #E: Healthy lung epithilium
]
vars = Dict(:labels=>[:E, ],
            :names=>["Healthy lung cells"],
            :initconds=>u0,
            :covars=>[1, ])
tspan = (0.0, 40.0)
pf = Dict(:labels=>[:a₁, :a₂, :b₁, :b₂, :α, :Vmin, :Vmax],
          :values=>[
            0.5     # a₁
            4.0     # a₂
            13.0    # b₁
            19.0    # b₂
            0.1     # α
            1e-8    # Vmin
            3.7e3   # Vmax
            ])
pp = Dict(:labels=>[:μ_E, :δ_E],
          :values=>[
            0.01    # μ_E
            0.15    # δ_E (x3) changed from notes to accommodate true V/Vmax and produce the same graph
            ])
pr = Dict(:labels=>[],
          :values=>[
            ])
pθ = Dict(:labels=>[:E, ],
          :values=>[
            0.3     # E_θ
            ])
pars = Dict(:fixed=>pf, :prior=>pp, :random=>pr, :thresholds=>pθ)
prob = ODEProblem(modelODE!, u0, tspan, parsDict2Array(pars));
sol = solve(prob, Tsit5(), dtmax=1e-1)
plot(sol)
model_VLF001 = RiskCov19ODEModel("VLF001", "VLF + Lung tissue", prob, vars, pars)

