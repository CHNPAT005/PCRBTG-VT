## Author: Patrick Chang
# Script file to simulate a multivariate Ornstein Uhlenbeck process

## Preamble

using Random, LinearAlgebra

function OU(n, μ, Σ, θ; kwargs...)
    # n - simulation length
    # mu - long term price average
    # sigma - covariance matrix of the stocks
    # theta - reversion parameter

    k = size(Σ)[1]

    kwargs = Dict(kwargs)

    if haskey(kwargs, :startprice)
        startprice = kwargs[:startprice]
    else
        startprice = fill(100.0, (k,1))
    end

    if haskey(kwargs, :seed)
        seed = kwargs[:seed]
    else
        seed = 1
    end

    P = zeros(n, 2)
    P[1,:] = log.(startprice)

    Random.seed!(seed)
    Z = randn(k, n-1)

    for i in 2:n
        P[i,:] = P[i-1,:] + θ .* (log.(μ) .- P[i-1,:]) + cholesky(sigma).L * Z[:,i-1]
    end

    return exp.(P)
end
