## Author: Patrick Chang
# Script file to investigate the correlation dynamic for 10 JSE equity
# using event time averaging. We consider skipping events that correspond to
# 1Hr, 10min and 1 min in Calendar time.

## Preamble

using CSV, Plots, Random, ProgressMeter, StatsBase, JLD, LaTeXStrings, DataFrames, Dates, Distributions

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
# Get the names of tickers
tickers = names(prices)

## Function to make event time data
#---------------------------------------------------------------------------

# Function to create the volume time data of with a specific number of buckets per day
function MakeEventTimeData(prices::DataFrame)
    # Remove all rows with NaNs
    inds = []
    for i in 1:size(prices)[1]
        if !all(isnan, prices[i,:])
            push!(inds, i)
        end
    end
    prices = prices[inds,:]
    # Get dimensions
    (n, m) = size(prices)
    # Create the new dataframe to store the bar data
    prices_ET = zeros(n, m)
    # Populate the first event. ASSUME: that events all start together, achieved by bringing the first obs in each ticker to the start
    for i in 1:m
        prices_ET[1,i] = filter!(!isnan, prices[:,i])[1]
    end
    # Loop through each line
    @showprogress "Computing..." for i in 2:n
        # Extract the row of data
        line = prices[i,:]
        for j in 1:m
            # Loop through each ticker
            if !isnan(line[j])
                # An event occured in that ticker => update
                prices_ET[i,j] = prices[i,j]
            else
                # No event occured in that ticker => keeps the current price
                prices_ET[i,j] = prices_ET[i-1,j]
            end
        end
    end
    return prices_ET, repeat(collect(1:1:n), 1, m)
end

# k-skip sampling such that a spefic number of events per day are retained
function kskip(data, numevents::Int64)
    # Get dimensions
    (n, m) = size(data)
    # Get skip size
    k = Int(floor(n / (5*numevents)))
    # Get the indicies
    inds = collect(1:k:n)
    # Extract the prices
    P = data[inds,:]
    return P, repeat(collect(1:k:n), 1, m)
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

## Full event data
EventData = MakeEventTimeData(prices)

MMED = NUFFTcorrDKFGG(convert(Matrix, EventData[1]), convert(Matrix, EventData[2]))[1]
HYED = RV(convert(Matrix, EventData[1]))[1]

## Event data, at a time scale ≡ to 1 min Calendar time
# Achieved using skip sampling for HY, but controlled using N in MM
SamplesPerday = 480
(n, m) = size(EventData[1])
Δt = Int(floor(n / (5*SamplesPerday)))

EventData480 = kskip(EventData[1], SamplesPerday)
N480 = Int(floor((n/Δt - 1)/2))

MMED480 = NUFFTcorrDKFGG(convert(Matrix, EventData[1]), convert(Matrix, EventData[2]), N = N480)[1]
HYED480 = RV(convert(Matrix, EventData480[1]))[1]

## Event data, at a time scale ≡ to 10 min Calendar time
# Achieved using skip sampling for HY, but controlled using N in MM
SamplesPerday = 48
(n, m) = size(EventData[1])
Δt = Int(floor(n / (5*SamplesPerday)))

EventData48 = kskip(EventData[1], SamplesPerday)
N48 = Int(floor((n/Δt - 1)/2))

MMED48 = NUFFTcorrDKFGG(convert(Matrix, EventData[1]), convert(Matrix, EventData[2]), N = N48)[1]
HYED48 = RV(convert(Matrix, EventData48[1]))[1]

## Event data, at a time scale ≡ to 1 Hr Calendar time
# Achieved using skip sampling for HY, but controlled using N in MM
SamplesPerday = 8
(n, m) = size(EventData[1])
Δt = Int(floor(n / (5*SamplesPerday)))

EventData8 = kskip(EventData[1], SamplesPerday)
N8 = Int(floor((n/Δt - 1)/2))

MMED8 = NUFFTcorrDKFGG(convert(Matrix, EventData[1]), convert(Matrix, EventData[2]), N = N8)[1]
HYED8 = RV(convert(Matrix, EventData8[1]))[1]

## Plot the results

# Function to compute the mean and standard deviation of the absolute value
# of the upper triangular matrix excluding the diagonal component
function summarise(ρ)
    D = size(ρ)[1]
    res = []
    for i in 1:(D-1)
        for j in (i+1):D
            res = [res; abs(ρ[i,j])]
        end
    end
    return mean(res), std(res)
end

## Event Time:
## 1 Hr equivalent
# MM
plot(MMED8, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :orange, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho^{ij}|} = $(round(summarise(MMED8)[1], digits = 4)) \\pm $(round(summarise(MMED8)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/MMED8.png")

# HY
plot(HYED8, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :orange, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho^{ij}|} = $(round(summarise(HYED8)[1], digits = 4)) \\pm $(round(summarise(HYED8)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/HYED8.png")


## 10 min equivalent
# MM
plot(MMED48, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :orange, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho^{ij}|} = $(round(summarise(MMED48)[1], digits = 4)) \\pm $(round(summarise(MMED48)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/MMED48.png")

# HY
plot(HYED48, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :orange, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho^{ij}|} = $(round(summarise(HYED48)[1], digits = 4)) \\pm $(round(summarise(HYED48)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/HYED48.png")


## 1 min equivalent
# MM
plot(MMED480, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :orange, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho^{ij}|} = $(round(summarise(MMED480)[1], digits = 4)) \\pm $(round(summarise(MMED480)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/MMED480.png")

# HY
plot(HYED480, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :orange, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho^{ij}|} = $(round(summarise(HYED480)[1], digits = 4)) \\pm $(round(summarise(HYED480)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/HYED480.png")

## TAQ equivalent
# MM
plot(MMED, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :orange, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho^{ij}|} = $(round(summarise(MMED)[1], digits = 4)) \\pm $(round(summarise(MMED)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/MMED.png")

# HY
plot(HYED, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :orange, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho^{ij}|} = $(round(summarise(HYED)[1], digits = 4)) \\pm $(round(summarise(HYED)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/HYED.png")
