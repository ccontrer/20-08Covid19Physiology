@model function turing_model(data::RiskCov19Data, ODEModel::RiskCov19ODEModel, ::Type{T} = Float64) where {T}
    # data: RiskCov19 structure containing (df, timecol, deadcol, covarcols, odegrcols)
    # ODEModel: RiskCov19 structure containing (id, name, prob, vars, pars)
    
    ## Info contained in data
    df = data.df
    N = size(df, 1) 
    t = df[:, data.timecol]
    d = df[:, data.deadcol]
    X = Matrix(df[:, data.lincovars])
    num_lin_covars = length(data.lincovars)
    Y = Matrix(df[:, data.odecovars])
    num_ode_groups = length(unique(Y))
            
    ## Info contained in ODEModel
    pars = ODEModel.pars
    vars = ODEModel.vars
    num_ode_covar = length(vars.covars_ix)
    num_rand_vars = length(pars.values) 

    ## Linear covariates
    μ_β = fill(0., num_lin_covars)
    Σ_β = diagm(vec(1 ./ std(X, dims=1)))
    β ~ MvNormal(μ_β, Σ_β)
    z = X*β
    
    ## ODE covariates
    # priors[i][k] is the value of the i-th parameter in the ODE of group k
    priors = Vector(undef, num_rand_vars)
    for i = 1:num_rand_vars
       priors[i] ~ filldist(truncated(Normal(pars.values[i], 0.01), 0.0, Inf), num_ode_groups)
    end
    η = Vector(undef, num_ode_groups)
    # sol[k] is the solution of the ODE for the k-th group
    sol = Vector(undef, num_ode_groups)
    # TODO: resolve t_s
    t_s = 5.0
    for k = 1:num_ode_groups
        η[k] = [priors[i][k] for i = 1:num_rand_vars]
        prob = remake(ODEModel.prob, p=η[k], tspan=(0.0, t_s + maximum(t)))
        sol[k] = solve(prob, Tsit5(), dt=0.1)
    end
    y(t, k) = sol[Int(k)+1](t)[vars.covars_ix]
    γ ~ filldist(Normal(0., 0.1), num_ode_covar)
    
    ## Survival model
    # Basal hazard rate
    # λ = 1.0
    λ ~ Exponential(1.0)
    if λ <= 0.0
        Turing.@addlogprob! -Inf
        return
    end
    # Exponential basal hazard function
    hazard₀(t, λ) = λ
    # hazard function
    hazard(t, z, k, λ) = hazard₀(t, λ)*exp(z + y(t, k)'*γ)
    
    # iterate for each patient
    for i in 1:N
        # cumulative hazard function
        h = hazard(t_s + t[i], z[i], Y[i], λ)
        H = try
            quadgk(s -> hazard(s, z[i], Y[i], λ), 0., t_s + t[i])[1]
        catch e
            Turing.@addlogprob! -Inf
            return
        end
        # TODO: add exception to unbounded integral
        # user Log-likelihood
        Turing.@addlogprob! d[i]*log(h) - H
    end
end

@model function turing_model(data::RiskCov19Data, ::Type{T} = Float64) where {T}
    # data: RiskCov19 structure containing (df, timecol, deadcol, covarcols, odegrcols)
    # ODEModel: RiskCov19 structure containing (id, name, prob, vars, pars)
    
    ## Info contained in data
    df = data.df
    N = size(df, 1) 
    t = df[:, data.timecol]
    d = df[:, data.deadcol]
    X = Matrix(df[:, data.lincovars])
    num_lin_covars = length(data.lincovars)
    Y = Matrix(df[:, data.odecovars])
    num_ode_groups = length(unique(Y))
            
    ## Info contained in ODEModel
    pars = ODEModel.pars
    vars = ODEModel.vars
    num_ode_covar = length(vars.covars_ix)
    num_rand_vars = length(pars.values) 

    ## Linear covariates
    μ_β = fill(0., num_lin_covars)
    Σ_β = diagm(vec(1 ./ std(X, dims=1)))
    β ~ MvNormal(μ_β, Σ_β)
    z = X*β
    
    ## Survival model
    # Basal hazard rate
    # λ = 1.0
    λ ~ Exponential(1.0)
    if λ <= 0.0
        Turing.@addlogprob! -Inf
        return
    end
    # Exponential basal hazard function
    hazard₀(t, λ) = λ
    # hazard function
    hazard(t, z, λ) = hazard₀(t, λ)*exp(z)
    
    # iterate for each patient
    for i in 1:N
        # cumulative hazard function
        h = hazard(t_s + t[i], z[i], λ)
        H = try
            quadgk(s -> hazard(s, z[i], λ), 0., t_s + t[i])[1]
        catch e
            Turing.@addlogprob! -Inf
            return
        end
        # TODO: add exception to unbounded integral
        # user Log-likelihood
        Turing.@addlogprob! d[i]*log(h) - H
    end
end

function rename_pars(chn::Chains, data::RiskCov19Data, model::RiskCov19ODEModel)
    num_ode_groups = length(unique(Matrix(data.df[:, data.odecovars])))
    names = String.(chn.name_map.parameters)
    filter!(x->occursin("priors", x), names)
    d = Dict{Symbol, Symbol}()
    for (i, name) in enumerate(String.(model.pars.names))
        for k = 1:num_ode_groups
            push!(d, Symbol("priors[$i][$k]") => Symbol("$name[$k]") )
        end
    end
    return replacenames(chn, d)
end

function train(
    data::RiskCov19Data, 
    model::RiskCov19ODEModel;
    n_samples=1_000, 
    n_chains=4,
    sampler=NUTS()
)
    m = turing_model(data, model)
    chn = mapreduce(c -> sample(m, sampler, n_samples), chainscat, 1:n_chains)
    return rename_pars(chn, data, model)
end

function train(
    data::RiskCov19Data;
    n_samples=1_000, 
    n_chains=4,
    sampler=NUTS()
)
    m = turing_model(data)
    chn = mapreduce(c -> sample(m, sampler, n_samples), chainscat, 1:n_chains)
    return rename_pars(chn, data, model)
end

# function solve_ODE_pred(chn::Chains, data::RiskCov19Data, model::RiskCov19ODEModel)
#     t = data.df[:, data.timecol]
#     num_ode_groups = length(unique(Matrix(data.df[:, data.odegrcols])))
#     par_names = string.(model.pars.names)
#     t_s = 5.0
#     sol = Vector{ODESolution}(undef, num_ode_groups)
#     for i = 1:num_ode_groups
#         p = [mean(chn, s) for s in Symbol.(par_names.*"[$i]")]
#         prob = remake(model.prob, p=p, tspan=(0.0, t_s + maximum(t)))
#         sol[i] = solve(prob, Tsit5(), dt=0.1)
#     end
#     return sol
# end

function predict(newdata::RiskCov19Data, chn::Chains, data::RiskCov19Data, model::RiskCov19ODEModel)
    
    ## Info contained in data
    df = newdata.df
    N = size(df, 1) 
    t = df[:, newdata.timecol]
    d = df[:, newdata.deadcol]
    X = Matrix(df[:, newdata.lincovars])
    Y = Matrix(df[:, newdata.odecovars])
    num_ode_groups = length(unique(Matrix(df[:, data.odecovars])))
    
    ## Info contained in ODEModel
    par_names = model.pars.names
    
    ## Estimated values
    # λ = mean(chn, :λ)
    λ = 0.1
    β = [mean(chn, s) for s in namesingroup(chn, :β)]
    γ = [mean(chn, s) for s in namesingroup(chn, :γ)]
    
    ## Linear covariates
    z = X*β
    
    ## Solve ODEs
    t_s = 5.0
    sol = Vector{ODESolution}(undef, num_ode_groups)
    for k = 1:num_ode_groups
        p = [mean(chn, s) for s in Symbol.(string.(par_names).*"[$k]")]
        prob = remake(model.prob, p=p, tspan=(0.0, t_s + maximum(t)))
        sol[k] = solve(prob, Tsit5(), dt=0.1)
    end
    y(t, k) = sol[Int(k)+1](t)[vars.covars_ix]
    
    ## Survival model
    # Exponential basal hazard function
    hazard₀(t, λ) = λ
    # hazard function
    hazard(t, z, k, λ) = hazard₀(t, λ)*exp(z + y(t, k)'*γ)
    
    ## Compute probability
    pred = Vector{Float64}(undef, N)
    for i = 1:N
        h = hazard(t_s + t[i], z[i], Y[i], λ)
        H = try
            quadgk(s -> hazard(s, z[i], Y[i], λ), 0., t_s + t[i])[1]
        catch e
            Inf
        end
        pred[i] = h*exp(-H)
    end
    return pred
end