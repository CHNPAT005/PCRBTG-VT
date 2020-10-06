## Author: Patrick Chang
# Script file to simulate multivariate Merton Model

## Preamble

using Random, LinearAlgebra, Distributions, PoissonRandom

## Building the Merton Model using the method of Paul Glasserman in his book
#   Monte Carlo Methods in Financial Engineering
#   This method allows me to build the path in 1 unit time increments
#   and there is no compensator in this method

function Merton(n, mu, sigma, lambda, a, b; kwargs...)
    # n - simulation length
    # mu - vector input of the drift component for the underlying GBM
    # sigma - covariance matrix of the stocks
    # lambda - vector input of the poisson process intensity
    # a - vector input of mean of normal jump size
    # b - vector input of std dev of normal jump size

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
    a = reshape(a, k, 1)
    b = reshape(b, k, 1)

    X = zeros(n, k)
    X[1,:] = log.(startprice)

    A = cholesky(sigma).L
    d = mu - sigma2./2 - lambda.*(exp.(a - 0.5*b.^2) .- 1)

    Random.seed!(seed)
    Z = randn(k, n-1)

    for i in 2:n
        M = zeros(k, 1)
        N = zeros(k, 1)

        for j in 1:k
            Random.seed!(seed+j+i+n)
            N[j] = pois_rand(lambda[j])
        end

        Random.seed!(seed+i)
        Z2 = randn(k)
        M = a.*N + b.*sqrt.(N).*Z2
        z = Z[:,i-1]
        X[i,:] = X[i-1,:] + d + A*z + M
    end
    return exp.(X)
end
