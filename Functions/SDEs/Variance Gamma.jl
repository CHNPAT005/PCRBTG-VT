## Author: Patrick Chang
# Script file to simulate a multivariate Variance Gamma

## Preamble

using Random, LinearAlgebra, Distributions

#---------------------------------------------------------------------------
## Building the Variance Gamma using the method of Paul Glasserman in his book
#   Monte Carlo Methods in Financial Engineering
#   This method allows me to build the path in 1 unit time increments

function VG(n, μ, Σ, β; kwargs...)
    # n - simlulation length
    # mu - vector input of the drift component
    # sigma - covariance matrix of the stocks
    # beta - scale for gamma

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

    P = zeros(n, k)
    P[1,:] = startprice

    Random.seed!(seed)
    Z = randn(k, n-1)
    Random.seed!(seed)
    Y = rand(Gamma(1/β, β), k, n-1)

    for i in 2:n
        z = Z[:,i-1]
        y = Y[:,i-1]
        P[i,:] = P[i-1,:] + μ .* y + cholesky(sigma).L * (sqrt.(y).*z)
    end

    return P
end
