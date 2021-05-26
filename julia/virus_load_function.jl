struct VLFResult
    # fit::Optim.MultivariateOptimizationResults
    fit::LsqFit.LsqFitResult
    param_array::AbstractArray
    data::VirusLoadData
    names::Vector
end

# TODO: implement confidence intervals

# H(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5)))
H(x) = 0.5*(tanh(1e3*x) + 1.)

v₁(t, a₁, a₂, logVmax) = 1. + (10^logVmax - 1.0)*(tanh(6.0*(t - (a₁ + a₂)/2)/abs(a₂ - a₁)) - tanh(-3.0*(a₂ + a₁)/abs(a₂ - a₁)))/2

v₂(t, a₂, α) = 1. - H(t - a₂) + H(t - a₂)*exp(-α*(t - a₂))

v₃(t, b₁, b₂, logVmin) = 1. - (1.0 - 10^logVmin)*(tanh(6.0*(t - (b₁ + b₂)/2)/abs(b₂ - b₁)) - tanh(-3.0*(b₂ + b₁)/abs(b₂ - b₁)))/2

function assert_params(a₁, a₂, b₁, b₂, α, logVmax)
    all([a₁ < a₂, a₂ < b₁, b₁ < b₂, a₁ > 0, a₂ > 0, b₁ > 0, b₂ > 0, α > 0, logVmax > 0])
end

function LogVirusLoadFunction(t, a₁, a₂, b₁, b₂, α, logVmax)
    logVmin = -7.0
    log10(v₁(t, a₁, a₂, logVmax)*v₂(t, a₂, α)*v₃(t, b₁, b₂, logVmin))
end

function LogEmpiricalVirusLoadFunction(t::Vector, p::Vector, data::VirusLoadData)
    θ = minimum(data.v)
    logV = LogVirusLoadFunction.(t, p...)
    max.(logV, θ)
end

function rss(p, data)
    sum(abs2, LogEmpiricalVirusLoadFunction(data.t, p, data) - data.v)
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
   
function get_rand_pars(fit1, data)
    tmax = maximum(data.t)
    p = fit1.param
    s = 1.1
    a₁ = rand(truncated(Normal(p[1], s), 0, tmax))
    a₂ = rand(truncated(Normal(p[2], s), a₁, tmax))
    b₁ = rand(truncated(Normal(p[3], s), a₂, tmax))
    b₂ = rand(truncated(Normal(p[4], s), b₁, tmax))
    α = rand(truncated(Normal(p[5], 0.5), 0, Inf))
    logVmax = rand(truncated(Normal(p[6], 0.5), 0, Inf))
    return [a₁, a₂, b₁, b₂, α, logVmax]
end

function fitVLF(data::VirusLoadData, p0::Vector)
    lb, ub = get_bounds(data, p0)
    return curve_fit((t, p) -> LogEmpiricalVirusLoadFunction(t, p, data), data.t, data.v, p0, lower=lb, upper=ub)
    # cost(p) = sum(abs2, LogEmpiricalVirusLoadFunction(data.t, p, data) - data.v)
    # fit = optimize(cost, lb, ub, p0, SAMIN(), Optim.Options(iterations=10^4))
end

function fitVLF(data::VirusLoadData; ϵ=0.1)
    # Best parameter
    fit1 = nothing
    niter = 10^3
    pb1 = Progress(niter, 0.5, "Fitting the VLF to data ")
    for iter = 1:niter
        #      a₁,     a₂,     b₁,     b₂      α,      logVmax
        par0 = vcat(sort(data.t[end]*rand(4)), rand(), maximum(data.v)*(rand()+0.5))
        fit0 = try 
            VirusLoadCurve.fitVLF(data, par0)
        catch y
            nothing
        end
        if fit0 != nothing
            if fit1 == nothing
                fit1 = fit0
            else
                if sum(abs2, fit0.resid) < sum(abs2, fit1.resid)
                    fit1 = fit0
                end
            end
        end
        next!(pb1)
    end
    # Set of parameters
    param_array = [fit1.param]
    niter = 10^6
    tol = (1.0 + ϵ)*rss(fit1.param, data)
    pb2 = Progress(niter, 0.5, "Finding possible parameter values ")
    for k = 1:niter
        par = get_rand_pars(fit1, data)
        if (rss(par, data) < tol && assert_params(par...))
            push!(param_array, par)
        end
        next!(pb2)
    end
    println("Number of possible parameter found: ", length(param_array))
    names = ["a₁", "a₂", "b₁", "b₂", "α", "logVmax"]
    VLFResult(fit1, param_array, data, names)
end

function param_extrema(param_array)
    if length(param_array)==0 return zip(fill(-Inf, 6), fill(Inf, 6)) end
    array = map(x -> map(y -> y[x], param_array), 1:length(param_array[1]))
    param_lower = map(minimum, array)
    param_upper = map(maximum, array)
    return zip(param_lower, param_upper)
end

@recipe function f(result::VLFResult; empirical=false, plotrange=true)
    tmin, tmax = extrema(result.data.t)
    tt = Vector(range(tmin, tmax, step=1e-3))
    func = empirical ? ((t, p) -> LogEmpiricalVirusLoadFunction(t, p, result.data)) : ((t, p) -> LogVirusLoadFunction.(t, p...))
    yy = func(tt, result.fit.param)
    vmin, vmax = extrema(result.data.v)
    xaxis := ("Time (days)", (tmin, tmax), font(14))
    yaxis := (L"\log\,V(t)", (vmin-0.5, vmax+0.5), font(14))
    grid --> :none
    
    @series begin
        x := tt
        y := yy
        linecolor := :black
        linewidth := 4
        label := ("Virus load function")
        ()
    end

    if plotrange
        yylower = yyupper = yy
        for param in result.param_array
            yy0 = func(tt, param)
            yylower = min.(yylower, yy0)
            yyupper = max.(yyupper, yy0)
        end
        vmax = max(vmax, maximum(yyupper))
        ylims := (vmin-0.5, vmax+0.5)
    
        @series begin
            x := tt
            y := yylower
            linewidth := 0
            linecolor := 2
            linealpha := 0.4
            fillcolor := 2
            fillalpha := 0.4
            fillrange := yyupper
            label := ("Range of VLF")
            ()
        end
    end

end

function Base.summary(result::VLFResult)
    println(@sprintf "RSS = %.5e" sum(abs2, result.fit.resid))
    [println(p, " = ", v, ", CI=", ci) for (p, v, ci) in zip(result.names, result.fit.param, param_extrema(result.param_array))]
    nothing
end
