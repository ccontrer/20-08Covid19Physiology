struct VTMResult
    fit::Any
    data::VirusLoadData
    names::Vector
    p0::Vector
    u0::Vector
end

ViralTargetODE = @ode_def begin
    dT = -β*T*V
    dI₁ = β*T*V - k*I₁
    dI₂ = k*I₁ - δ*I₂/(K + I₂)
    dV = p*I₂ - c*V
end β k δ K p c

u0default = [1e+7, 75.0, 0.0, 1e-12]
lbdefault = [0.0, 0.0, 0.0, 0.0, 0.0]
ubdefault = [Inf, Inf, Inf, Inf, Inf]

function LogViralTargetModel(t, p, data, u0)
    tspan = (0.0, maximum(data.t))
    pars = (β = p[1],
            k = 4.0,
            δ = p[2],
            K = p[3],
            p = p[4],
            c = p[5])
    prob = ODEProblem(ViralTargetODE, u0, tspan, pars)
    sol = solve(prob, Tsit5(), dtmax=1e-2)
    log10.(sol(t)[end, :])
end

function LogEmpiricalViralTargetModel(t, p, data, u0)
    θ = minimum(data.v)
    LogV = LogViralTargetModel(t, p, data, u0)
    max.(LogV, θ)
end

function fitVTM(data::VirusLoadData, p0::Vector; lb=lbdefault, ub=ubdefault, u0=u0default)
    model(t, p) = LogEmpiricalViralTargetModel(t, p, data, u0)
    fit = curve_fit(model, data.t, data.v, p0, lower=lb, upper=ub)
    names = ["β", "δ", "K", "p", "c"]
    VTMResult(fit, data, names, p0, u0)
end

@recipe function f(result::VTMResult; empirical=false, p0=false)
    tmin, tmax = 0.0, maximum(result.data.t)
    vmin, vmax = extrema(result.data.v)
    param = p0 ? result.fit.param : result.p0
    tt = Vector(range(tmin, tmax, step=1e-2))
    x := tt
    if empirical
        y := LogEmpiricalViralTargetModel(tt, param, result.data, result.u0)
    else
        y := LogViralTargetModel(tt, param, result.data, result.u0)
    end
    linewidth --> 4
    linestyle --> :dash
    xaxis --> ("Time (days)", (tmin, tmax))
    yaxis --> (L"\log\,V(t)", (vmin-1, vmax+1))
    grid --> :none
    label --> "Viral target model"
    ()
end

function Base.summary(result::VTMResult)
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