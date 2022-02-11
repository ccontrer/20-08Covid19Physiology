

# basal hazard function
# TODO: generalize parametric distribution
# h₀(t, dist) = pdf(dist, t)/(ccdf(dist, t))
h₀(t, dist) = params(dist)[1]

# hazard function
h(t, z, dist) = h₀(t, dist)*exp(z)

# cumulative hazard function
# TODO: optimize integration
H(t, z, dist) = quadgk(s -> h(s, z, dist), 0., t)[1]

# log-likelihood function
# TODO: resolve dist datatype and ForwardDifff
logL(t::T, d::T, z, dist::ContinuousUnivariateDistribution) where {T<:Real} = 
    d*log(h(t, z, dist)) - H(t, z, dist)