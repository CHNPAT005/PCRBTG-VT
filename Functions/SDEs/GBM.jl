## Author: Patrick Chang
# Script file to simulate a multivariate Geometric Brownian Motion

## Preamble

using Random, LinearAlgebra

## Code for Geometric Brownian Motion using the method from Paul Glasserman
#  in his book - Monte Carlo Methods in Financial Engineering

# timesteps (dt) can be controlled by scaling mu and sigma accordingly

function GBM(n, mu, sigma; kwargs...)
    # n - simlulation length
    # mu - vector input of the drift component
    # sigma - covariance matrix of the stocks
    # startprice - starting price for the assets

    # all inputs must have appropriate dimensions

    k = size(sigma)[1]

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

    mu = reshape(mu, k, 1)
    sigma = reshape(sigma, k, k)
    sigma2 = reshape(diag(sigma), k, 1)

    P = zeros(n, k)
    P[1,:] = startprice

    A = cholesky(sigma).L
    b = mu - sigma2./2

    Random.seed!(seed)
    Z = randn(k, n-1)

    for i in 2:n
        z = Z[:,i-1]
        X = b + A * z
        P[i,:] = P[i-1,:] .* exp.(X)
    end
    return P
end
