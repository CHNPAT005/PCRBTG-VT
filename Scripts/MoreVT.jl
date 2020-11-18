## Author: Patrick Chang
# Script file to further investigate the Epps effect under volume time
# averaging, with different distributional assumptions for volume

## Preamble

using Pipe, StatsBase, ProgressMeter, Distributions, JLD, Distributions

cd("/Users/patrickchang1/PCRBTG-VT")

# Hawkes
include("../Functions/Hawkes/Hawkes.jl")

#---------------------------------------------------------------------------
## Supporting functions to streamline the process

# Epps stuff
#------------------------
# Theoretical limiting correlations from the Hawkes price model used in Barcy et al.
function theoreticalCorr(α_12, α_13, β)
    Γ_12 = α_12 / β; Γ_13 = α_13 / β
    num = 2*Γ_13*(1+Γ_12)
    den = 1 + Γ_13^2 + 2*Γ_12 + Γ_12^2

    return num/den
end
# Theoretical Epps effect in Calendar time from the Hawkes price model used in Barcy et al.
function theoreticalEpps(τ, μ, α_12, α_13, β)
    Γ_12 = α_12 / β; Γ_13 = α_13 / β
    Λ = μ / (1 - Γ_12 - Γ_13)
    Q_1 = -(μ * (Γ_12^2 + Γ_12 - Γ_13^2)) / (((Γ_12 + 1)^2 - Γ_13^2) * (1 - Γ_12 - Γ_13))
    Q_2 = -(μ * Γ_13) / (((Γ_12 + 1)^2 - Γ_13^2) * (1 - Γ_12 - Γ_13))
    R = (β * μ) / (Γ_12 + Γ_13 - 1)
    G_1 = β * (1 + Γ_12 + Γ_13)
    G_2 = β * (1 + Γ_12 - Γ_13)
    C_1 = (2 + Γ_12 + Γ_13) * (Γ_12 + Γ_13) / (1 + Γ_12 + Γ_13)
    C_2 = (2 + Γ_12 - Γ_13) * (Γ_12 - Γ_13) / (1 + Γ_12 - Γ_13)

    C_11 = Λ + (R*C_1)/(2*G_1) + (R*C_2)/(2*G_2) + R * (C_2*G_1^2*exp(-τ*G_2) - C_1*G_2^2 + Q_1*G_2^2*exp(-τ*G_1) -C_2*G_1^2) / (2*G_2^2*G_1^2*τ)
    C_12 = - (R*C_1)/(2*G_1) + (R*C_2)/(2*G_2) + R * (C_1*G_2^2 - C_2*G_1^2 - C_1*G_2^2*exp(-τ*G_1) + C_2*G_1^2*exp(-τ*G_2)) / (2*G_2^2*G_1^2*τ)

    return C_12/C_11
end

# Extracting the prices
#------------------------
# Construct Prices from the simulation paths
function makePrices(P0, uptick, downtick)
    alltimes = [uptick; downtick]
    alltimes = sort(alltimes)

    P = []
    for i in 1:length(alltimes)
        up = findall(uptick .<= alltimes[i])
        down = findall(downtick .<= alltimes[i])
        p = P0 + length(up) - length(down)
        P = append!(P, p)
    end
    return P, alltimes
end

# Helpful functions
#------------------------
# Realised Volatility estimate
# Used in place of HY when observations are synchronous (as they are the same), for computational efficiency
function RV(prices)
    δ = diff(log.(prices), dims = 1)
    Σ = δ' * δ
    σ = sqrt.(diag(Σ))
    ρ = Σ ./ (σ * σ')
    return ρ, Σ
end
# Simulate Pareto distribution
function RPareto(xn, α, n)
    return xn ./ ((rand(n)).^(1/α))
end

# Making the prices
#------------------------
# Construct volume time price path with uniform distribution
function VolumeTimeUniformPrice(Asset1, Asset2, num::Int, maxvol::Int, seed)
    ### num: the number of samples we want
    # Get prices from each asset
    P1 = exp.(Asset1[1]); P2 = exp.(Asset2[1])
    # Simulate random volumes to a max of maxvol
    Random.seed!(seed)
    V1 = sample(1:maxvol, length(P1), replace = true)
    Random.seed!(seed + length(P1))
    V2 = sample(1:maxvol, length(P2), replace = true)
    # Expand price obs by volume sample
    P1 = vcat(fill.(P1, Int.(V1))...)
    P2 = vcat(fill.(P2, Int.(V2))...)
    # Get size of volume buckets
    V1size = Int(floor(length(P1)/num))
    V2size = Int(floor(length(P2)/num))
    # Initialize storage
    P = zeros(num,2)
    # Loop through the data
    for i in 1:num
        P[i,1] = mean(P1[(i-1)*V1size+1:i*V1size])
        P[i,2] = mean(P2[(i-1)*V2size+1:i*V2size])
    end
    return log.(P)
end
# Construct volume time price path with normal distribution
function VolumeTimeGaussPrice(Asset1, Asset2, num::Int, seed)
    ### num: the number of samples we want
    # Get prices from each asset
    P1 = exp.(Asset1[1]); P2 = exp.(Asset2[1])
    # Simulate random volumes to a max of maxvol
    Random.seed!(seed)
    V1 = round.(rand(Normal(50, 5), length(P1)))
    Random.seed!(seed + length(P1))
    V2 = round.(rand(Normal(50, 5), length(P2)))
    # Expand price obs by volume sample
    P1 = vcat(fill.(P1, Int.(V1))...)
    P2 = vcat(fill.(P2, Int.(V2))...)
    # Get size of volume buckets
    V1size = Int(floor(length(P1)/num))
    V2size = Int(floor(length(P2)/num))
    # Initialize storage
    P = zeros(num,2)
    # Loop through the data
    for i in 1:num
        P[i,1] = mean(P1[(i-1)*V1size+1:i*V1size])
        P[i,2] = mean(P2[(i-1)*V2size+1:i*V2size])
    end
    return log.(P)
end
# Construct volume time price path with pareto distribution
function VolumeTimeParetoPrice(Asset1, Asset2, num::Int, seed)
    ### num: the number of samples we want
    # Get prices from each asset
    P1 = exp.(Asset1[1]); P2 = exp.(Asset2[1])
    # Simulate random volumes to a max of maxvol
    Random.seed!(seed)
    V1 = round.(RPareto(10, 1.7, length(P1)))
    Random.seed!(seed + length(P1))
    V2 = round.(RPareto(10, 1.7, length(P2)))
    # Expand price obs by volume sample
    P1 = vcat(fill.(P1, Int.(V1))...)
    P2 = vcat(fill.(P2, Int.(V2))...)
    # Get size of volume buckets
    V1size = Int(floor(length(P1)/num))
    V2size = Int(floor(length(P2)/num))
    # Initialize storage
    P = zeros(num,2)
    # Loop through the data
    for i in 1:num
        P[i,1] = mean(P1[(i-1)*V1size+1:i*V1size])
        P[i,2] = mean(P2[(i-1)*V2size+1:i*V2size])
    end
    return log.(P)
end
# Construct volume time price path with beta distribution
function VolumeTimeBetaPrice(Asset1, Asset2, num::Int, scale, seed)
    ### num: the number of samples we want
    # Get prices from each asset
    P1 = exp.(Asset1[1]); P2 = exp.(Asset2[1])
    # Simulate random volumes to a max of maxvol
    Random.seed!(seed)
    V1 = ceil.(rand(Beta(scale), length(P1)) .* 100)
    Random.seed!(seed + length(P1))
    V2 = ceil.(rand(Beta(scale), length(P2)) .* 100)
    # Expand price obs by volume sample
    P1 = vcat(fill.(P1, Int.(V1))...)
    P2 = vcat(fill.(P2, Int.(V2))...)
    # Get size of volume buckets
    V1size = Int(floor(length(P1)/num))
    V2size = Int(floor(length(P2)/num))
    # Initialize storage
    P = zeros(num,2)
    # Loop through the data
    for i in 1:num
        P[i,1] = mean(P1[(i-1)*V1size+1:i*V1size])
        P[i,2] = mean(P2[(i-1)*V2size+1:i*V2size])
    end
    return log.(P)
end

#---------------------------------------------------------------------------
## Main functions to obtain the results

# Function to extract estimates under volume time over different bucket sizes
function getVolumeTimeCorrelations(Asset1, Asset2, dt, T, seed)
    # Get range of time scales
    n = length(dt)
    # Initialize storage variables for the estimators
    Unif500 = zeros(n, 1)
    Unif200 = zeros(n, 1)
    Unif100 = zeros(n, 1)
    Unif50 = zeros(n, 1)
    Gauss = zeros(n, 1)
    Pareto = zeros(n, 1)
    B1 = zeros(n, 1)
    B2 = zeros(n, 1)
    B10 = zeros(n, 1)
    # Loop through the ranges
    for i in 1:n
        # Create volume data for given certain amount of samples
        Uni500 = VolumeTimeUniformPrice(Asset1, Asset2, Int(floor(T/dt[i])), 500, seed)
        Uni200 = VolumeTimeUniformPrice(Asset1, Asset2, Int(floor(T/dt[i])), 200, seed)
        Uni100 = VolumeTimeUniformPrice(Asset1, Asset2, Int(floor(T/dt[i])), 100, seed)
        Uni50 = VolumeTimeUniformPrice(Asset1, Asset2, Int(floor(T/dt[i])), 50, seed)
        Gau = VolumeTimeGaussPrice(Asset1, Asset2, Int(floor(T/dt[i])), seed)
        Par = VolumeTimeParetoPrice(Asset1, Asset2, Int(floor(T/dt[i])), seed)
        Be1 = VolumeTimeBetaPrice(Asset1, Asset2, Int(floor(T/dt[i])), 0.1, seed)
        Be2 = VolumeTimeBetaPrice(Asset1, Asset2, Int(floor(T/dt[i])), 0.2, seed)
        Be10 = VolumeTimeBetaPrice(Asset1, Asset2, Int(floor(T/dt[i])), 1, seed)
        # Compute the results
        Unif500[i] = RV(exp.(Uni500))[1][1,2]
        Unif200[i] = RV(exp.(Uni200))[1][1,2]
        Unif100[i] = RV(exp.(Uni100))[1][1,2]
        Unif50[i] = RV(exp.(Uni50))[1][1,2]
        Gauss[i] = RV(exp.(Gau))[1][1,2]
        Pareto[i] = RV(exp.(Par))[1][1,2]
        B1[i] = RV(exp.(Be1))[1][1,2]
        B2[i] = RV(exp.(Be2))[1][1,2]
        B10[i] = RV(exp.(Be10))[1][1,2]
    end
    return Unif500, Unif200, Unif100, Unif50, Gauss, Pareto, B1, B2, B10
end

# Function to strealine everything together
function getRes(reps; maxvol = 200)
    # Set up parameters
    T = 3600*20
    par_1 = BarcyParams(0.015, 0.023, 0.05, 0.11)
    lambda0_1 = par_1[1]; alpha_1 = par_1[2]; beta_1 = par_1[3]
    domain = collect(1:1:100)
    Random.seed!(2020)
    seeds = Int.(floor.(rand(reps) .* 1000000))
    # Initialize storage variables
    U500 = zeros(reps, length(domain)); U200 = zeros(reps, length(domain)); U100 = zeros(reps, length(domain))
    U50 = zeros(reps, length(domain)); Gauss = zeros(reps, length(domain)); Pareto = zeros(reps, length(domain))
    B1 = zeros(reps, length(domain)); B2 = zeros(reps, length(domain)); B10 = zeros(reps, length(domain))
    # Loop through the different cases
    @showprogress "Computing..." for i in 1:reps
        # Simulate the processs
        t = simulateHawkes(lambda0_1, alpha_1, beta_1, T, seed = seeds[i])
        # Obtain the results
        VolumeTimeRes = getVolumeTimeCorrelations(makePrices(0, t[1], t[2]), makePrices(0, t[3], t[4]), domain, T, seeds[i])
        # Store the results
        U500[i,:] = VolumeTimeRes[1]; U200[i,:] = VolumeTimeRes[2]; U100[i,:] = VolumeTimeRes[3]
        U50[i,:] = VolumeTimeRes[4]; Gauss[i,:] = VolumeTimeRes[5]; Pareto[i,:] = VolumeTimeRes[6]
        B1[i,:] = VolumeTimeRes[7]; B2[i,:] = VolumeTimeRes[8]; B10[i,:] = VolumeTimeRes[9]
    end
    return U500, U200, U100, U50, Gauss, Pareto, B1, B2, B10
end

## Obtain the results
reps = 100
res = getRes(reps)

# Save results
save("Computed Data/MoreVT.jld", "res", res)
# Load results
res = load("Computed Data/MoreVT.jld")
res = res["res"]

## Plot the results
dt = collect(1:1:100)
ρ = theoreticalCorr(0.023, 0.05, 0.11)
q = quantile.(TDist(reps-1), [0.975])

# Volume time (all)
p1 = plot(dt, mean(res[1], dims=1)', ribbon = (q .* std(res[1], dims=1)'), fillalpha = .1, legend = :topleft, line=(1, [:dash]), label = L"\textrm{U500}", dpi = 300)
plot!(p1, dt, mean(res[2], dims=1)', ribbon = (q .* std(res[2], dims=1)'), fillalpha=.1, line=(1, [:dash]), label = L"\textrm{U200}")
plot!(p1, dt, mean(res[3], dims=1)', ribbon = (q .* std(res[3], dims=1)'), fillalpha=.1, line=(1, [:dash]), label = L"\textrm{U100}")
plot!(p1, dt, mean(res[4], dims=1)', ribbon = (q .* std(res[4], dims=1)'), fillalpha=.1, line=(1, [:dash]), label = L"\textrm{U50}")
plot!(p1, dt, mean(res[5], dims=1)', ribbon = (q .* std(res[5], dims=1)'), fillalpha=.1, line=(1, [:dash]), label = L"\textrm{Gauss}")
plot!(p1, dt, mean(res[6], dims=1)', ribbon = (q .* std(res[6], dims=1)'), fillalpha=.1, line=(1, [:dash]), label = L"\textrm{Pareto}")
plot!(p1, dt, mean(res[7], dims=1)', ribbon = (q .* std(res[7], dims=1)'), fillalpha=.1, line=(1, [:dash]), label = L"\textrm{B.1}")
plot!(p1, dt, mean(res[8], dims=1)', ribbon = (q .* std(res[8], dims=1)'), fillalpha=.1, line=(1, [:dash]), label = L"\textrm{B.2}")
plot!(p1, dt, mean(res[9], dims=1)', ribbon = (q .* std(res[9], dims=1)'), fillalpha=.1, line=(1, [:dash]), label = L"\textrm{B1}")
xlabel!(p1, L"\textrm{Sampling interval (volume time)}")
ylabel!(p1, L"\tilde{\rho}_{\Delta t}^{ij}")


p2 = plot(dt, mean(res[3], dims=1)', ribbon = (q .* std(res[1], dims=1)'), fillalpha = .1, legend = :topleft, line=(1, [:dash]), label = L"\textrm{Uniform(1, 100)}")
plot!(p2, dt, mean(res[5], dims=1)', ribbon = (q .* std(res[2], dims=1)'), fillalpha=.1, line=(1, [:dash]), label = L"\textrm{Normal(50, 5)}")
plot!(p2, dt, mean(res[6], dims=1)', ribbon = (q .* std(res[3], dims=1)'), fillalpha=.1, line=(1, [:dash]), label = L"\textrm{Power-law}")
plot!(p2, dt, mean(res[7], dims=1)', ribbon = (q .* std(res[4], dims=1)'), fillalpha=.1, line=(1, [:dash]), label = L"\textrm{Beta(0.1, 0.1)}")
plot!(p2, dt, mean(res[8], dims=1)', ribbon = (q .* std(res[5], dims=1)'), fillalpha=.1, line=(1, [:dash]), label = L"\textrm{Beta(0.2, 0.2)}")
plot!(p2, dt, mean(res[9], dims=1)', ribbon = (q .* std(res[6], dims=1)'), fillalpha=.1, line=(1, [:dash]), label = L"\textrm{Beta(2, 2)}")
xlabel!(p2, L"\textrm{Sampling interval (volume time)}")
ylabel!(p2, L"\tilde{\rho}_{\Delta t}^{ij}")
savefig(p2, "Plots/MoreVT.svg")
