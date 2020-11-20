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

function LogViralTargetModel(t, p, data, u0)
    θ = minimum(data.v)
    tspan = (0.0, maximum(data.t))
    pars = (β = p[1],
            k = p[2],
            δ = p[3],
            K = p[4],
            p = p[5],
            c = p[6])
    prob = ODEProblem(ViralTargetODE, u0, tspan, pars)
    sol = solve(prob, Tsit5(), dtmax=1e-2)
    max.(log10.(sol(t)[end, :]), θ)
end

function fitVTM(data::VirusLoadData, p0::Vector; u0=u0default)
    model(t, p) = LogViralTargetModel(t, p, data, u0)
    lb = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    ub = [Inf, Inf, Inf, Inf, Inf, Inf]
    fit = curve_fit(model, data.t, data.v, p0, lower=lb, upper=ub)
    names = ["β", "k", "δ", "K", "p", "c"]
    VTMResult(fit, data, names, p0, u0)
end

@recipe function f(result::VTMResult)
    tmin, tmax = 0.0, maximum(result.data.t)
    vmin, vmax = extrema(result.data.v)
    tt = Vector(range(tmin, tmax, step=1e-2))
    x := tt
    y := LogViralTargetModel(tt, result.fit.param, result.data, result.u0)
    linewidth := 4
    linestyle := :dash
    xaxis --> ("Time (days)", (tmin, tmax))
    yaxis --> (L"\log\,V(t)", (vmin-1, vmax+1))
    grid --> :none
    label := "Viral target model"
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