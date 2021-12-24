### Normal Tissue Comlication Model ODE models
### v0.1

struct NTCM_model_ODE
    id::String
    name::String
    prob::ODEProblem
    u_labels::NamedTuple
    u_cofactor_ix::Vector
    pars::Dict
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
u_labels = [:E, ]
u_names = ["Healthy lung cells"]
u_ntuple = (; zip(u_labels, u_names)...)
u0 = [
    1.0   #E: Healthy lung epithilium
]
tspan = (0.0, 40.0)
pf_labels = [:a₁, :a₂, :b₁, :b₂, :α, :Vmin, :δ_E];
pf = [
    0.5     # a₁
    4.0     # a₂
    13.0    # b₁
    19.0    # b₂
    0.1     # α
    1e-8    # Vmin
    0.15    # δ_E (x3) changed from notes to accommodate true V/Vmax and produce the same graph
    ]
pf_ntuple = (; zip(pf_labels, pf)...)
pp_labels = [];
pp = [
    ]
pp_ntuple = (; zip(pp_labels, pp)...)
pr_labels = [:Vmax, :μ_E];
pr = [
    3.7e3   # Vmax
    0.01    # μ_E
    ]
pr_ntuple = (; zip(pr_labels, pr)...)
pθ_labels = [:E, ];
pθ = [
    0.3     # E_θ
    ]
pθ_ntuple = (; zip(pθ_labels, pθ)...)
pars = Dict(:fixed=>pf_ntuple, :prior=>pp_ntuple, :random=>pr_ntuple, :thresholds=>pθ_ntuple)
prob = ODEProblem(modelODE!, u0, tspan, vcat(pf, pp, pr));
# prob = ODEProblem(modelODE!, u0, tspan, vcat(collect(values(pars[:fixed])), collect(values(pars[:prior])), collect(values(pars[:rand]))));
sol = solve(prob, Tsit5(), dtmax=1e-1)
plot(sol)
model_ODE_VLF001 = NTCM_model_ODE("VLF001", "VLF + Lung tissue", prob, u_ntuple, pars)


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
