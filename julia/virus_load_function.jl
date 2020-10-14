struct VLFResult
    fit::Any
    data::VirusLoadData
    names::Vector
    p0::Vector
end

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5)))

v₁(t, a₁, a₂, logVmax) = 1. + (10^logVmax - 1.0)*(tanh(6.0*(t - (a₁ + a₂)/2)/(a₂ - a₁)) - tanh(-3.0*(a₂ + a₁)/(a₂ - a₁)))/2

v₂(t, a₂, α) = 1. - heaviside(t - a₂) + heaviside(t - a₂)*exp(-α*(t - a₂))

v₃(t, b₁, b₂, logVmin) = 1. - (1.0 - 10^logVmin)*(tanh(6.0*(t - (b₁ + b₂)/2)/(b₂ - b₁)) - tanh(-3.0*(b₂ + b₁)/(b₂ - b₁)))/2

function assert_params(a₁, a₂, b₁, b₂, α, logVmax)
    if a₁ > a₂ @error("a₁ > a₂") end
    if a₂ > b₁ @error("a₂ > b₁") end
    if b₁ > b₂ @error("b₁ > b₂") end
end

function LogVirusLoadFunction(t, p, data::VirusLoadData)
    θ = minimum(data.v)
    logVmin = -6.0
    a₁, a₂, b₁, b₂, α, logVmax = p
    logV = log10.(v₁.(t, a₁, a₂, logVmax).*v₂.(t, a₂, α).*v₃.(t, b₁, b₂, logVmin))
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
    model(t, p) = LogVirusLoadFunction(t, p, data)
    lb, ub = get_bounds(data, p0)
    fit = curve_fit(model, data.t, data.v, p0, lower=lb, upper=ub)
    names = ["a₁", "a₂", "b₁", "b₂", "α", "logVmax"]
    VLFResult(fit, data, names, p0)
end

@recipe function f(result::VLFResult)
    tmin, tmax = 0.0, maximum(result.data.t)
    vmin, vmax = extrema(result.data.v)
    tt = Vector(range(tmin, tmax, step=1e-2))
    x := tt
    y := LogVirusLoadFunction(tt, result.fit.param, result.data)
    linewidth := 3
    xaxis --> ("Time (days)", (tmin, tmax))
    yaxis --> (L"\log\,V(t)", (vmin-1, vmax+1))
    grid --> :none
    label := "Virus load function"
    ()
end

function Base.summary(result::VLFResult)
    println(@sprintf "RSS = %.5e (convergence: %s)" sum(abs2, result.fit.resid) result.fit.converged)
    CI = try 
        confidence_interval(result.fit, 0.05)
    catch
        false
    end
    for i in 1:length(result.fit.param)
        name = result.names[i]
        val = result.fit.param[i]
        p0 = result.p0[i]
        if CI == false
            println(@sprintf "  %s = %.3e (initial=%.3e)" name val p0)
        
        else
            CIl = CI[i][1]
            CIr = CI[i][2]
            println(@sprintf "  %s = %.3e (CI=(%.3e, %.3e), initial=%.3e)" name val CIl CIr p0)
        end
    end
end
