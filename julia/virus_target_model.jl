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

du[1]  = dAngI = PRA - (c_ACE + c_NEP)*AngI - ln(2.)/h_AngI*AngI
du[2]  = dAngII = c_ACE*AngI - k_ACE2*ACE2m*AngII - (c_AT1R + c_AT2R)*AngII - ln(2.)*AngII/h_AngII
du[3]  = dAT1AngII = c_AT1R*AngII - ln(2.)*AT1RAngII/h_AT1R
du[4]  = dAT2AngII = c_AT2R*AngII - ln(2.)*AT2RAngII/h_AT2R
du[5]  = dAng17 = c_NEP*AngI + k_ACE2*ACE2m*AngII - ln(2.)*Ang17/h_Ang17
du[6]  = dACE2m = S_ACE2m - (c_ADAM17 + c_AT1R*AT1RAngII + c_Cov19*P_plasma/(K_mACE2i + P_plasma + c_ACE2i*ACE2p))*ACE2m - ln(2.)*ACE2m/h_ACE2
du[7]  = dACE2p = c_ADAM17*ACE2m - ln(2.)*ACE2p/h_ACE2
du[8]  = dACE2i = (c_AT1R*AT1RAngII + c_Cov19*P_plasma/(K_mACE2i + P_plasma + c_ACE2i*ACE2p))*ACE2m - ln(2.)*ACE2i/h_ACE2
du[9]  = dP_plasma = k_burst*P - k_ACE2P*ACE2p*P_plasma     - c_Cov19*P_plasma*ACE2m/(K_mACE2i + P_plasma _ c_ACE2i*ACE2p)
du[10] = dP = k_Cov19*c_Cov19*P_plasma*ACE2m/(K_mACE2i + P_plasma _ c_ACE2i*ACE2p) + k_pg*P*(1 - P/P_∞) - k_pm*s_m*P/(μ_m + k_mp*P) - k_pn*f(N) - k_burst*P
du[11] = dN = s_nr*R/(μ_nr + R) + k_AT1R*AT1RAngII - μ_n*N
du[12] = dD = k_dn*f_s(f(N)) - μ_d*D
du[13] = dC_A = s_c + k_cn*Q/(1 + Q) + k_ACE2p*ACE2p + k_Ang17*Ang17 + k_AT2R*AT2RAngII - μ_c*C_A