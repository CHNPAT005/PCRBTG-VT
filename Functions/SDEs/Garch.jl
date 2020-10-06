## Author: Patrick Chang
# Script file to simulate a bivariate GARCH (1,1)

## Preamble

using Random, LinearAlgebra

#---------------------------------------------------------------------------
## Building the bivariate GARCH (1,1) as specified by
# Roberto Reno - 2003

function Garch(n, θ, λ, ω, ρ; kwargs...)
    # Kwarg setup
    kwargs = Dict(kwargs)

    if haskey(kwargs, :startprice)
        startprice = kwargs[:startprice]
    else
        startprice = log.(fill(100.0, (2,1)))
    end

    if haskey(kwargs, :seed)
        seed = kwargs[:seed]
    else
        seed = 1
    end

    if haskey(kwargs, :startvar)
        startvar = kwargs[:startvar]
    else
        startvar = [0.5; 0.6]
    end

    if haskey(kwargs, :dt)
        dt = kwargs[:dt]
    else
        dt = 86400
    end

    # Extract the variables
    λ1 = λ[1]; λ2 = λ[2]
    ω1 = ω[1]; ω2 = ω[2]
    θ1 = θ[1]; θ2 = θ[2]

    Random.seed!(seed)
    Z = randn(4, n-1)

    ## Sigma2 path
    Sigma2 = zeros(2, n)
    Sigma2[:,1] = startvar

    for i in 2:n
        Sigma2[1,i] = Sigma2[1,i-1] + λ1 * (ω1 - Sigma2[1,i-1])/dt + sqrt(2*λ1*θ1*Sigma2[1,i-1]/dt)*Z[3, i-1]
        Sigma2[2,i] = Sigma2[2,i-1] + λ2 * (ω2 - Sigma2[2,i-1])/dt + sqrt(2*λ2*θ2*Sigma2[2,i-1]/dt)*Z[4, i-1]
    end

    P = zeros(n, 2)
    P[1,:] = startprice

    for i in 2:n
        Σ = [Sigma2[1,i] sqrt(Sigma2[1,i] * Sigma2[2,i])*ρ; sqrt(Sigma2[1,i] * Sigma2[2,i])*ρ Sigma2[2,i]]
        P[i,:] = cholesky(sigma).L * Z[1:2, i-1] / sqrt(dt) + P[i-1,:]
    end

    return exp.(P)
end

theta = [0.035, 0.054]
lambda = [0.296, 0.48]
w = [0.636, 0.476]

P = Garch(n, theta, lambda, w, 0.35)

NUFFTcorrDKFGG(P, t)[1]
