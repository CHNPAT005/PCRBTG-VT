## Author: Patrick Chang
# Script file to investigate the Epps effect in volume time
# for 4 banking stocks on the JSE.

## Preamble

using CSV, Plots, Random, ProgressMeter, StatsBase, JLD, LaTeXStrings, DataFrames, Dates, Distributions, Pipe, ColorSchemes

cd("/Users/patrickchang1/PCRBTG-VT")

# Estimators
include("../Functions/Correlation Estimators/Dirichlet/NUFFTcorrDK-FGG.jl")
include("../Functions/Correlation Estimators/HY/HYcorr.jl")

# Read in the data
prices = CSV.read("Real Data/JSE_prices_2019-06-24_2019-06-28.csv")
times = CSV.read("Real Data/JSE_times_2019-06-24_2019-06-28.csv")
volume = CSV.read("Real Data/JSE_volume_2019-06-24_2019-06-28.csv")

# Remove extra column
prices = prices[:,2:end]; times = times[:,2:end]; volume = volume[:,2:end]
# Pull out banking stocks
prices = prices[:,[6:8; 10]]; times = times[:,[6:8; 10]]; volume = volume[:,[6:8; 10]]
# Get the names of tickers
tickers = names(prices)

#---------------------------------------------------------------------------
## Supporting functions to streamline the process
# Realised Volatility estimate
# Used in place of HY when observations are synchronous (as they are the same), for computational efficiency
function RV(prices)
    δ = diff(log.(prices), dims = 1)
    Σ = δ' * δ
    σ = sqrt.(diag(Σ))
    ρ = Σ ./ (σ * σ')
    return ρ, Σ
end
# Function to split the data into the 5 days
# Results in dictionary of data, with a vector (of the two assets) of matricies
# where the matrix is price, times, volume
function DataSplit(A1::String, A2::String, P::DataFrame, t::DataFrame, V::DataFrame)
    # Filter out the pair of interest
    p1 = filter(!isnan, P[:,A1]); p2 = filter(!isnan, P[:,A2])
    t1 = filter(!isnan, t[:,A1]); t2 = filter(!isnan, t[:,A2])
    V1 = filter(!isnan, V[:,A1]); V2 = filter(!isnan, V[:,A2])
    # Convert the times to dates for index extraction later on
    A1dates = Date.(unix2datetime.(t1))
    A2dates = Date.(unix2datetime.(t2))
    dates_unique = unique(A1dates)
    # Initialize storage
    data = Dict()
    # Loop through each day
    for i in 1:length(dates_unique)
        # Extract the data
        date_indsA1 = findall(x -> x == dates_unique[i], A1dates)
        date_indsA2 = findall(x -> x == dates_unique[i], A2dates)
        # Data for the day
        day_data = []
        push!(day_data, [p1[date_indsA1] t1[date_indsA1] V1[date_indsA1]])
        push!(day_data, [p2[date_indsA2] t2[date_indsA2] V2[date_indsA2]])
        # Add to dictionary
        push!(data, i => day_data)
    end
    return data
end
# Function to make the volume time data into RV format
function MakeRVData(data, Δt; T = 28200)
    # Extract prices and volume
    P1 = data[1][:,1]; P2 = data[2][:,1]
    V1 = data[1][:,3]; V2 = data[2][:,3]
    # Expand price obs by volume sample
    P1 = vcat(fill.(P1, Int.(V1))...)
    P2 = vcat(fill.(P2, Int.(V2))...)
    # Get number of samples required
    num = Int(floor(T/Δt))
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
    return P
end
# Function to get the volume time correlations
function getVTcorrs(tickers; P=prices, t=times, V=volume)
    # Compute the number of pairwise comparisons
    npairs = Int(factorial(length(tickers)) / (factorial(2) * factorial(length(tickers)-2)))
    # Time scale of investigation
    dt = collect(1:1:300)
    # Initialize the estimates
    RVest = zeros(length(dt), npairs)
    # Set up ind for storage
    ind = 1
    # Loop through the pairs
    @showprogress "Computing..." for i in 1:(length(tickers)-1)
        for j in (i+1):length(tickers)
            # Split the data into separate days
            data = DataSplit(tickers[i], tickers[j], P, t, V)
            # Initialize temporary storage matricies for the estimates
            RVtemp = zeros(length(dt), 5)
            # Loop through the different days
            for k in 1:length(data)
                # Extract data for the day
                day_data = data[k]
                # Loop through the different time scales
                for l in 1:length(dt)
                    # Make RV data
                    RVData = MakeRVData(day_data, dt[l])
                    # Compute correlations
                    RVtemp[l,k] = RV(RVData)[1][1,2]
                end
            end
            # Store the 5 day average into the matrix
            RVest[:,ind] = mean(RVtemp, dims = 2)
            # Update the ind
            ind += 1
        end
    end
    return RVest
end

## Obtain the results
res = getVTcorrs(tickers)

# Save results
save("Computed Data/EmpVT.jld", "res", res)
# Load results
res = load("Computed Data/EmpVT.jld")
res = res["res"]

## Extract label pairs
pairnames = Matrix{Union{Nothing, String}}(nothing, 1, 6)
inds = 1
for i in 1:(length(tickers)-1)
    for j in (i+1):length(tickers)
        # push!(pairnames, "$(tickers[i])"*"/"*"$(tickers[j])")
        pairnames[inds] = "$(tickers[i])"*"/"*"$(tickers[j])"
        inds += 1
    end
end

## Plot the results
dt = collect(1:1:300)

plot(dt, res, legend = :topleft, palette = ColorSchemes.tab10.colors, label = pairnames)
xlabel!(L"\textrm{Sampling interval (volume time)}")
ylabel!(L"\tilde{\rho}_{\Delta t}^{ij}")
savefig("Plots/EmpVT.svg")

## Obtain the results for appendix
# Read in the data
prices = CSV.read("Real Data/JSE_prices_2019-06-24_2019-06-28.csv")
times = CSV.read("Real Data/JSE_times_2019-06-24_2019-06-28.csv")
volume = CSV.read("Real Data/JSE_volume_2019-06-24_2019-06-28.csv")

# Remove extra column
prices = prices[:,2:end]; times = times[:,2:end]; volume = volume[:,2:end]
# Pull out banking stocks
prices = prices[:,[1;3;6;7]]; times = times[:,[1;3;6;7]]; volume = volume[:,[1;3;6;7]]
# Get the names of tickers
tickers = names(prices)

res = getVTcorrs(tickers)

# Save results
save("Computed Data/EmpAppVT.jld", "res", res)
