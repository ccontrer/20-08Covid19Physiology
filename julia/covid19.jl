using OrdinaryDiffEq, ParameterizedFunctions
using Optim
using Plots, LaTeXStrings
using DelimitedFiles
using Statistics
using Revise

import LsqFit
include("julia/NLFit.jl")
import .NLFit

tdata = [range(0, 12 , step=1.0);]
tend = tdata[end]
data = readdlm("../data/Smith2018/Virus_Best10.txt", ' ', Float64, '\n')'

ttdata = reshape(repeat(tdata', 10), 1, :)[:]
vvdata = reshape(data', 1, :)[:]

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5)))
v₁(t, a₁, a₂, max) = 1. + (max - 1.0)*(tanh(6.0*(t - (a₁ + a₂)/2)/(a₂ - a₁)) - tanh(-3.0*(a₂ + a₁)/(a₂ - a₁)))/2
v₂(t, a₂, α) = 1. - heaviside(t - a₂) + heaviside(t - a₂)*exp(-α*(t - a₂))
v₃(t, b₁, b₂, min) = 1. - (1.0 - min)*(tanh(6.0*(t - (b₁ + b₂)/2)/(b₂ - b₁)) - tanh(-3.0*(b₂ + b₁)/(b₂ - b₁)))/2
function ViralLoad(t, par) 
    a₁, a₂, b₁, b₂, α, mini, maxi = par
#     log10.(v₁.(t, a₁, a₂, maxi).*v₂.(t, a₂, α).*v₃.(t, b₁, b₂, mini))
    # TODO: fix negative numbers
    log10.(max.(v₁.(t, a₁, a₂, maxi).*v₂.(t, a₂, α).*v₃.(t, b₁, b₂, mini), 1e-15))
end

par0 = [1.0, 5.3, 6.0, 8.5, 0.5, 1e-6, 1.5e6]
lower = [0.0, 0.0, 0.0, 0.0, 1e-10, 1e-12, 1e2]
upper = [12.0, 12.0, 12.0, 12., 1.0, 1e2, 1e7]
inner_optimizer = ConjugateGradient()
opts = Optim.Options(g_tol=1e-2)
fit = NLFit.curve_fit(ViralLoad, ttdata, vvdata, lower, upper, par0, Fminbox(inner_optimizer), opts)

fit2 = LsqFit.curve_fit(ViralLoad, ttdata, vvdata, par0)



println("Estimated parameter values:")
println("   \tRSS       \ta_1 \ta_2 \tb_1 \tb_2 \talpha \tmin \tmax")
println(@sprintf "a \t%f" 1.0)
print(@sprintf " Init.\t%4.2e\t%2.3f% \t%2.3f \t%2.3f \t%2.3f \t%2.3f \t%1.0e \t%1.0e " sum(abs2, fit.resid) fit.parsam)
# print(" Est. \t{:4.2e}\t{:2.3f} \t{:2.3f} \t{:2.3f} \t{:2.3f} \t{:2.3f} \t{:1.0e} \t{:1.0e} ".format(
#     vl_avg.RSS, *vl_avg.par))
# print(" SE     \t     \t%1.0e \t%1.0e \t%1.0e \t%1.0e \t%1.0e \t%1.0e \t%1.0e" % tuple(vl_avg.par_se))

# https://github.com/JuliaNLSolvers/Optim.jl/blob/master/docs/src/examples/maxlikenlm.jl