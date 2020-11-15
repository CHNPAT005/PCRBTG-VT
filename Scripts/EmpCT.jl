## Author: Patrick Chang
# Script file to investigate the Epps effect in calendar time
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
# Function to make the calendar time data into MM and HY format
# pad the data and pull the prices and times together
function MakeMMHYData(data)
    # Pull out the data
    A1 = data[1]; A2 = data[2]
    # Get the dimensions
    n1 = size(A1, 1); n2 = size(A2, 1)
    n = max(n1,n2)
    # Initialize price and time matrix
    P = fill(NaN, n, 2); t = fill(NaN, n, 2)
    # Populate the price and time matrix
    P[1:n1,1] = A1[:,1]; P[1:n2,2] = A2[:,1]
    t[1:n1,1] = A1[:,2]; t[1:n2,2] = A2[:,2]
    return P, t
end
# Function to make the calendar time data into RV format, apply previous tick interpolation
function MakeRVData(data, Δt; T = 28200)
    # Pull out the data
    A1 = data[1]; A2 = data[2]
    # Extract the prices and time
    p1 = A1[:,1]; p2 = A2[:,1]
    t1 = A1[:,2]; t2 = A2[:,2]
    # Re-scale the times
    T0 = min(t1[1], t2[1])
    t1 = t1 .- T0; t2 = t2 .- T0
    # For the lagging asset, bring the first observation backwards so that both assets start at the same time
    if t1[1] != 0
        t1 = [0; t1]
        p1 = [p1[1]; p1]
    end
    if t2[1] != 0
        t2 = [0; t2]
        p2 = [p2[1]; p2]
    end
    # Sample times
    st = collect(0:Δt:T)
    # Initialize the PTI price matrix
    n = length(st)
    PTI = zeros(n, 2)
    # Loop through all the sample times
    for i in 1:length(st)
        # PTI indicies
        ind1 = findlast(x -> x <= st[i], t1)
        ind2 = findlast(x -> x <= st[i], t2)
        # Store into PTI matrix
        PTI[i,1] = p1[ind1]
        PTI[i,2] = p2[ind2]
    end
    return PTI
end
# Function to get the calendar time correlations
function getCTcorrs(tickers; P=prices, t=times, V=volume)
    # Compute the number of pairwise comparisons
    npairs = Int(factorial(length(tickers)) / (factorial(2) * factorial(length(tickers)-2)))
    # Time scale of investigation
    dt = collect(1:1:300)
    # Initialize the estimates
    MM = zeros(length(dt), npairs)
    RVest = zeros(length(dt), npairs)
    HY = zeros(npairs, 1)
    # Set up ind for storage
    ind = 1
    # Loop through the pairs
    @showprogress "Computing..." for i in 1:(length(tickers)-1)
        for j in (i+1):length(tickers)
            # Split the data into separate days
            data = DataSplit(tickers[i], tickers[j], P, t, V)
            # Initialize temporary storage matricies for the estimates
            MMtemp = zeros(length(dt), 5)
            RVtemp = zeros(length(dt), 5)
            HYtemp = zeros(5, 1)
            # Loop through the different days
            for k in 1:length(data)
                # Extract data for the day
                day_data = data[k]
                # Create the MM and HY dataset
                MMHYData = MakeMMHYData(day_data)
                # Loop through the different time scales
                for l in 1:length(dt)
                    # Make RV data
                    RVData = MakeRVData(day_data, dt[l])
                    # Get N for MM
                    N = Int(floor((28200/dt[l] - 1)/2))
                    # Compute correlations
                    MMtemp[l,k] = NUFFTcorrDKFGG(MMHYData[1], MMHYData[2], N = N)[1][1,2]
                    RVtemp[l,k] = RV(RVData)[1][1,2]
                end
                # Compute HY corr
                HYtemp[k] = HYcorr(MMHYData[1][:,1], MMHYData[1][:,2], MMHYData[2][:,1], MMHYData[2][:,2])[1][1,2]
            end
            # Store the 5 day average into the matrix
            MM[:,ind] = mean(MMtemp, dims = 2)
            RVest[:,ind] = mean(RVtemp, dims = 2)
            HY[ind]   = mean(HYtemp)
            # Update the ind
            ind += 1
        end
    end
    return MM, RVest, HY
end


## Obtain the results
res = getCTcorrs(tickers)

# Save results
save("Computed Data/EmpCT.jld", "res", res)
# Load results
res = load("Computed Data/EmpCT.jld")
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

plot(dt, res[1], legend = :bottomright, palette = ColorSchemes.tab10.colors, label = pairnames, ylims = (0,0.82))
xlabel!(L"\textrm{Sampling interval (calendar time)}")
ylabel!(L"\tilde{\rho}_{\Delta t}^{ij}")
savefig("Plots/EmpCTMM.svg")


plot(dt, res[2], legend = :bottomright, palette = ColorSchemes.tab10.colors, label = pairnames, ylims = (0,0.82))
xlabel!(L"\textrm{Sampling interval (calendar time)}")
ylabel!(L"\tilde{\rho}_{\Delta t}^{ij}")
savefig("Plots/EmpCTRV.svg")


plot(dt, repeat([NaN], 300), legend = :topright, ylims = (0,0.65), label = "")
hline!(res[3]', palette = ColorSchemes.tab10.colors, label = pairnames)
xlabel!(L"\textrm{Sampling interval (calendar time)}")
ylabel!(L"\tilde{\rho}_{\Delta t}^{ij}")
savefig("Plots/EmpCTHY.svg")
