### Normal Tissue Comlication Model ODE models
### v0.1

struct NTCM_model
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
    H(t) = 0.5*(tanh(1e3*t) + 1.)
    v₁(t, a₁, a₂, Vmax) = 1. + (Vmax - 1.0)*(tanh(6.0*(t - (a₁ + a₂)/2)/(a₂ - a₁)) - tanh(-3.0*(a₂ + a₁)/(a₂ - a₁)))/2
    v₂(t, a₂, α) = 1 - H(t - a₂) + H(t - a₂)*exp(-α*(t - a₂))
    v₃(t, b₁, b₂, Vmin) = 1 - (1.0 - Vmin)*(tanh(6.0*(t - (b₁ + b₂)/2)/(b₂ - b₁)) - tanh(-3.0*(b₂ + b₁)/(b₂ - b₁)))/2
    v₁(t, a₁, a₂, Vmax)*v₂(t, a₂, α)*v₃(t, b₁, b₂, Vmin)
end

## Model: VLF001
## VLF + Lung tissue
# ODE model
function modelODE!(du, u, p, t)
    # p = [pf, pp, pr]
    E, = u
    a₁, a₂, b₁, b₂, α, Vmin, δ_E, Vmax, μ_E = p
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
pf = Dict(:labels=>[:a₁, :a₂, :b₁, :b₂, :α, :Vmin, :δ_E],
          :values=>[
            0.5     # a₁
            4.0     # a₂
            13.0    # b₁
            19.0    # b₂
            0.1     # α
            1e-8    # Vmin
            0.15    # δ_E (x3) changed from notes to accommodate true V/Vmax and produce the same graph
            ])
pp = Dict(:labels=>[],
          :values=>[
            ])
pr = Dict(:labels=>[:Vmax, :μ_E],
          :values=>[
            3.7e3   # Vmax
            0.01    # μ_E
            ])
pθ = Dict(:labels=>[:E, ],
          :values=>[
            0.3     # E_θ
            ])
pars = Dict(:fixed=>pf, :prior=>pp, :random=>pr, :thresholds=>pθ)
prob = ODEProblem(modelODE!, u0, tspan, parsDict2Array(pars));
sol = solve(prob, Tsit5(), dtmax=1e-1)
plot(sol)
model_VLF001 = NTCM_model("VLF001", "VLF + Lung tissue", prob, vars, pars)


# ## VLF + Lung + Cytokine + Immune
# # ODE model
# function modelODE!(du, u, p, t)
#     E, C, X₁, X₂ = u
#     a₁, a₂, b₁, b₂, α, Vmax, Vmin, μ_E, δ_E, μ_C, M_C, γ_C, α_X, γ_X, β₁, β₂, E_θ = p
#     V(t) = VLF(t, [a₁, a₂, b₁, b₂, α, Vmax, Vmin])
#     # Lung damage
#     du[1] = dE = μ_E*E*(1.0 - E) - δ_E*E*V(t)/Vmax # Healthy lung epithelium
#     du[2] = dC = μ_C*(M_C - C) + γ_C*V(t)*(1 - E)
#     du[3] = dX₁ = α_X*C*X₁ - γ_X*V(t)/Vmax*X₁ - β₁*X₁
#     du[4] = dX₂ = γ_X*V(t)/Vmax*X₁ - β₂*X₂
# end
# u_labels = [:E, :C, :X₁, :X₂]
# u_names = ["Healthy lung cells" "Cytokines" "Healthy mmune cells" "Infected immune cells"]
# u_ntuple = (; zip(u_labels, u_names)...)
# u0 = [
#     1.0   #E: Healthy lung epithilium
#     0.001 #C: Cytokines promoters
#     2.0   #X₁: Immune cells
#     0.0   #X₂: Infected immune cells
# ]
# tspan = (0.0, 40.0)
# p_labels = [:a₁, :a₂, :b₁, :b₂, :α, :Vmax, :Vmin, :μ_E, :δ_E, :μ_C, :M_C, :γ_C, :α_X, :γ_X, :β₁, :β₂, :E_θ];
# p = [
#     0.5     # a₁
#     4.0     # a₂
#     13.0    # b₁
#     19.0    # b₂
#     0.1     # α
#     3.7e3    # Vmax
#     1e-8    # Vmin
#     0.01    # μ_E 
#     0.15    # δ_E (x3) changed from notes to accommodate true V/Vmax and produce the same graph
#     0.5     # μ_C
#     0.1     # M_C
#     0.01    # γ_C
#     0.01    # α_X
#     0.000001   # γ_X (x1e-3) changed from notes to accommodate true V/Vmax and produce the same graph
#     0.00001    # β₁ (x1e-3) changed from notes to accommodate true V/Vmax and produce the same graph
#     0.5     # β₂
#     0.3     # E_θ
#     ]
# p_ntuple = (; zip(p_labels, p)...)
# prob = ODEProblem(modelODE!, u0, tspan, p);
