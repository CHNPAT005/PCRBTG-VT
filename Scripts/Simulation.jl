## Author: Patrick Chang
# Script file to investigate the Epps effect in event time
# using a Hawkes process to simulate the prices.
# We use the price model by Bacry et al. (2013) and compare
# the MM and HY estimator.

## Preamble

using Pipe, StatsBase, ProgressMeter, Distributions

cd("/Users/patrickchang1/PCRBTG-VT")

# Estimators
include("../Functions/Correlation Estimators/Dirichlet/NUFFTcorrDK-FGG.jl")
include("../Functions/Correlation Estimators/HY/HYcorr.jl")

# Hawkes
include("../Functions/Hawkes/Hawkes.jl")

#---------------------------------------------------------------------------
## Supporting functions to streamline the process
# Theoretical correlations from a Hawkes price model used in Barcy et al.
function theoreticalCorr(α_12, α_13, β)
    Γ_12 = α_12 / β; Γ_13 = α_13 / β
    num = 2*Γ_13*(1+Γ_12)
    den = 1 + Γ_13^2 + 2*Γ_12 + Γ_12^2

    return num/den
end

# Construct Prices from the simulation paths
function makePrices(P0, uptick, downtick)
    alltimes = [uptick; downtick]
    alltimes = sort(alltimes)

    P = []
    for i in 1:length(alltimes)
        up = findall(uptick .<= alltimes[i])
        down = findall(downtick .<= alltimes[i])
        p = P0 + 0.1 *(length(up) - length(down))
        P = append!(P, p)
    end
    return P, alltimes
end
# Construct event time price path
function EventTimePrice(Asset1, Asset2, P01, P02)
    # Get trading times and prices from each asset
    t1 = Asset1[2]; t2 = Asset2[2]
    P1 = Asset1[1]; P2 = Asset2[1]
    # Combine all trading times
    times = [0;t1;t2]
    sort!(times)
    # Initialize variables
    n = length(times)
    p1 = zeros(n)
    p2 = zeros(n)
    # Starting value
    p1[1] = P01
    p2[1] = P02
    # Counters
    i_1 = 1
    i_2 = 1
    # Loop through each event
    for i in 2:n
        if in(times[i], t1)
            p1[i] = P1[i_1]
            p2[i] = p2[i-1]
            i_1 += 1
        else
            p2[i] = P2[i_2]
            p1[i] = p1[i-1]
            i_2 += 1
        end
    end
    return [p1 p2], [collect(0.0:1:n-1) collect(0.0:1:n-1)]
end
# Construct volume time price path
function VolumeTimePrice(prices, num::Int)
    # Get number of price observations
    n = size(prices)[1]
    # Initialize volume averaged prices
    vbucket = Int(floor(n/num))
    # nn = @pipe collect(1:num:n) |> length
    P = zeros(vbucket, 2)
    # Loop through the data
    for i in 1:vbucket
        inds = (i-1)*num+1:i*num
        P[i,:] = @pipe inds |> prices[_,:] |> mean(_, dims = 1)
    end
    return P
end
# k-skip sampling
function kskip(prices, k)
    # Get number of price observations
    n = size(prices)[1]
    # Get the indicies
    inds = collect(1:k:n)
    # Extract the prices
    P = prices[inds,:]
    return P
end
# Realised Volatility estimate
# Used in place of HY when observations are synchronous (as they are the same)
# Used for computational efficiency
function RV(prices)
    δ = diff(log.(prices), dims = 1)
    Σ = δ' * δ
    σ = sqrt.(diag(Σ))
    ρ = Σ ./ (σ * σ')
    return ρ, Σ
end

#---------------------------------------------------------------------------
## Main functions to obtain the results

# Function to extract estimates under event time over different sampling intervals
function getEventTimeCorrelations(data, dt)
    # Get range of time scales
    prices = data[1]
    t = data[2]
    n = length(dt)
    # Initialize storage variables for the estimators
    MM = zeros(n, 1)
    HY = zeros(n, 1)
    # Loop through the ranges
    for i in 1:n
        # MM time scale
        N = Int(floor((size(prices)[1]/dt[i] - 1)/2))
        # HY k-skip
        tempdata = kskip(prices, dt[i])
        # Compute the results
        MM[i] = NUFFTcorrDKFGG(exp.(prices), t, N = N)[1][1,2]
        HY[i] = RV(exp.(tempdata))[1][1,2]
    end
    return MM, HY
end

# Function to extract estimates under volume time over different bucket sizes
function getVolumeTimeCorrelations(data, bucketsizes)
    # Get range of time scales
    prices = data[1]
    # t = data[2]
    n = length(bucketsizes)
    # Initialize storage variables for the estimators
    MM = zeros(n, 1)
    HY = zeros(n, 1)
    # Loop through the ranges
    for i in 1:n
        # HY k-skip
        tempdata = VolumeTimePrice(prices, bucketsizes[i])
        t = [collect(1:1:size(tempdata)[1]) collect(1:1:size(tempdata)[1])]
        # Compute the results
        MM[i] = NUFFTcorrDKFGG(exp.(tempdata), t)[1][1,2]
        HY[i] = RV(exp.(tempdata))[1][1,2]
    end
    return MM, HY
end

# Function to strealine everything together
function getRes(reps)
    # Set up parameters
    T = 3600*20
    par_1 = BarcyParams(0.015, 0.023, 0.05, 0.11)
    lambda0_1 = par_1[1]; alpha_1 = par_1[2]; beta_1 = par_1[3]
    domain = collect(1:1:100)
    Random.seed!(2020)
    seeds = Int.(floor.(rand(reps) .* 1000000))
    # Initialize storage variables
    MMET = zeros(reps, length(domain)); MMVT = zeros(reps, length(domain))
    HYET = zeros(reps, length(domain)); HYVT = zeros(reps, length(domain))
    # Loop through the different cases
    @showprogress "Computing..." for i in 1:reps
        # Simulate the processs
        t = simulateHawkes(lambda0_1, alpha_1, beta_1, T, seed = seeds[i])
        # Get raw event times
        data = EventTimePrice(makePrices(0, t[1], t[2]), makePrices(0, t[3], t[4]), 0, 0)
        # Obtain the results
        EventTimeRes = getEventTimeCorrelations(data, domain)
        VolumeTimeRes = getVolumeTimeCorrelations(data, domain)
        # Store the results
        MMET[i,:] = EventTimeRes[1]; HYET[i,:] = EventTimeRes[2]
        MMVT[i,:] = VolumeTimeRes[1]; HYVT[i,:] = VolumeTimeRes[2]
    end
    return MMET, MMVT, HYET, HYVT
end


## Obtain the results
reps = 100
test = getRes(reps)


## Plot the results
ρ = theoreticalCorr(0.023, 0.05, 0.11)
q = quantile.(TDist(reps-1), [0.975])


p1 = plot(1:100, mean(test[1], dims=1)', ribbon = (q .* std(test[1], dims=1)'), fillalpha = .3, legend = :bottomright, color = :blue, line=(1, [:solid]), label = L"\textrm{MM}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)), dpi = 300)
plot!(p1, dt, mean(test[3], dims=1)', ribbon = (q .* std(test[3], dims=1)'), fillalpha=.3, color = :red, line=(1, [:dash]), label = L"\textrm{HY}", marker = ([:x :d], 3, 0.8, stroke(3, :red)))
hline!(p1, [ρ], color = :black, line=(2, [:dot]), label = L"\textrm{Limiting } \rho")
xlabel!(p1, L"\textrm{Event time } \Delta t")
ylabel!(p1, L"\tilde{\rho}_{\Delta t}^{ij}")
# savefig("Plots/SimET.svg")


p2 = plot(1:100, mean(test[2], dims=1)', ribbon = (q .* std(test[2], dims=1)'), fillalpha = .3, legend = :bottomright, color = :blue, line=(1, [:solid]), label = L"\textrm{MM}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)), dpi = 300)
plot!(p2, dt, mean(test[4], dims=1)', ribbon = (q .* std(test[4], dims=1)'), fillalpha=.3, color = :red, line=(1, [:dash]), label = L"\textrm{HY}", marker = ([:x :d], 3, 0.8, stroke(3, :red)))
hline!(p2, [ρ], color = :black, line=(2, [:dot]), label = L"\textrm{Limiting } \rho")
xlabel!(p2, L"\textrm{Size of volume bucket } \Delta v")
ylabel!(p2, L"\tilde{\rho}_{\Delta v}^{ij}")
# savefig("Plots/SimVT.svg")


test = EventTimePrice(makePrices(0, t[1], t[2]), makePrices(0, t[3], t[4]), 0, 0)
plot(test[2][1:10,:], test[1][1:10,:], seriestype = :steppost)
