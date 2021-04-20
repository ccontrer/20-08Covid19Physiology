struct VLFResult
    fit::Optim.MultivariateOptimizationResults
    data::VirusLoadData
    names::Vector
    p0::Vector
end

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5)))

v₁(t, a₁, a₂, logVmax) = 1. + (10^logVmax - 1.0)*(tanh(6.0*(t - (a₁ + a₂)/2)/abs(a₂ - a₁)) - tanh(-3.0*(a₂ + a₁)/abs(a₂ - a₁)))/2

v₂(t, a₂, α) = 1. - heaviside(t - a₂) + heaviside(t - a₂)*exp(-α*(t - a₂))

v₃(t, b₁, b₂, logVmin) = 1. - (1.0 - 10^logVmin)*(tanh(6.0*(t - (b₁ + b₂)/2)/abs(b₂ - b₁)) - tanh(-3.0*(b₂ + b₁)/abs(b₂ - b₁)))/2

function assert_params(a₁, a₂, b₁, b₂, α, logVmax)
    if a₁ > a₂ @error("a₁ > a₂") end
    if a₂ > b₁ @error("a₂ > b₁") end
    if b₁ > b₂ @error("b₁ > b₂") end
end

function LogVirusLoadFunction(t, p, data::VirusLoadData)
    logVmin = -7.0
    a₁, a₂, b₁, b₂, α, logVmax = p
    log10.(v₁.(t, a₁, a₂, logVmax).*v₂.(t, a₂, α).*v₃.(t, b₁, b₂, logVmin))
end

function LogEmpiricalVirusLoadFunction(t, p, data::VirusLoadData)
    θ = minimum(data.v)
    logV = LogVirusLoadFunction(t, p, data::VirusLoadData)
    max.(logV, θ)
end

function get_bounds(data, p0)
    a₁, a₂, b₁, b₂, α, logVmax = p0
    a₁a₂ = (a₁ + a₂)/2
    a₂b₁ = (a₂ + b₁)/2
    b₁b₂ = (b₁ + b₂)/2
    tmin = 0.0
    tmax = maximum(data.t)
    vmax = maximum(data.v)
    lb   = [tmin, a₁a₂, a₂b₁, b₁b₂, 1e-8, vmax-3.0]
    ub   = [a₁a₂, a₂b₁, b₁b₂, tmax, 1e+4, vmax+3.0]
    return lb, ub
end
    

function fitVLF(data::VirusLoadData, p0::Vector)
    cost(p) = sum(abs2, LogEmpiricalVirusLoadFunction(data.t, p, data) - data.v)
    lb, ub = get_bounds(data, p0)
    fit = optimize(cost, lb, ub, p0, SAMIN(), Optim.Options(iterations=10^4))
    names = ["a₁", "a₂", "b₁", "b₂", "α", "logVmax"]
    VLFResult(fit, data, names, p0)
end

@recipe function f(result::VLFResult; empirical=false, stderrors=false)
    tmin, tmax = 0.0, maximum(result.data.t)
    vmin, vmax = extrema(result.data.v)
    tt = Vector(range(tmin, tmax, step=1e-2))
    x := tt
    if empirical
        yy = LogEmpiricalVirusLoadFunction(tt, result.fit.minimizer, result.data)
    else
        yy = LogVirusLoadFunction(tt, result.fit.minimizer, result.data)
    end
    y := yy
    # if stderrors
    #     CI = try 
    #         confidence_interval(result.fit, 0.05)
    #     catch
    #         false
    #     end
    #     yymin = yymax = yy
    #     if ~(CI == false)
    #         for k in 1:50
    #             par0 = [(CI[i][2]-CI[i][1])*rand()+CI[i][1] for i in 1:length(CI)]
    #             try 
    #                 if empirical
    #                     yy0 = LogEmpiricalVirusLoadFunction(tt, par0, result.data)
    #                 else
    #                     yy0 = LogVirusLoadFunction(tt, par0, result.data)
    #                 end
    #                 yymin = min.(yymin, yy0)
    #                 yymax = max.(yymax, yy0)
    #             catch
    #                 continue
    #             end
    #         end
    #     end
    #     ribbon := (yy-yymin, yymax-yy)
    # end
    linewidth --> 4
    xaxis --> ("Time (days)", (tmin, tmax))
    yaxis --> (L"\log\,V(t)", (vmin-1, vmax+1))
    grid --> :none
    label --> "Virus load function"
    ()
end

function Base.summary(result::VLFResult)
    println(@sprintf "RSS = %.5e" result.fit.minimum)
    # CI = try 
    #     confidence_interval(result.fit, 0.05)
    # catch
    #     false
    # end
    for i in 1:length(result.fit.minimizer)
        name = result.names[i]
        val = result.fit.minimizer[i]
        p0 = result.p0[i]
        # if CI == false
        println(@sprintf "  %s = %.3e (initial=%.3e)" name val p0)
        # else
        #     CIl = CI[i][1]
        #     CIr = CI[i][2]
        #     println(@sprintf "  %s = %.3e (CI=(%.3e, %.3e), initial=%.3e)" name val CIl CIr p0)
        # end
    end
end
