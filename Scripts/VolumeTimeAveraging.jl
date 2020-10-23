## Author: Patrick Chang
# Script file to investigate the correlation dynamic for 10 JSE equity
# using volume time averaging. We consider bucket sizes that correspond to
# 1Hr, 10min and 1 min in Calendar time. There are bucket sizes of
# 8, 48 and 480 buckets per day.

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

## Function to make volume time data
#---------------------------------------------------------------------------

# Function to create the volume time data of with a specific number of buckets per day
function MakeVolumeTimeData(prices::DataFrame, volume::DataFrame, numbuckets::Integer)
    # Create the new dataframe to store the bar data
    prices_VT = zeros(numbuckets*5, size(prices)[2])
    # Loop through each column
    @showprogress "Computing..." for i in 1:size(prices)[2]
        # Extract volume and prices for a particular asset
        temp_price = filter!(!isnan, prices[:,i])
        temp_vol = filter!(!isnan, volume[:,i])
        # Repeat each price as many times as the number of volume
        long_prices = vcat(fill.(temp_price, Int.(temp_vol))...)
        # Get the buket sizes
        bucketsize = Int(floor(sum(temp_vol)/(numbuckets*5)))
        # Loop through all the buckets
        for j in 1:(numbuckets*5)
            prices_VT[j,i] = mean(long_prices[(j-1)*bucketsize+1:j*bucketsize])
        end
    end
    return prices_VT, repeat(collect(1:1:numbuckets*5), 1, size(prices)[2])
end

## 8 buckets per day ≡ 1 Hour in Calendar time
VB8 = MakeVolumeTimeData(prices, volume, 8)

MMVB8 = NUFFTcorrDKFGG(convert(Matrix, VB8[1]), convert(Matrix, VB8[2]))[1]
HYVB8 = HYcorrFull(convert(Matrix, VB8[1]), convert(Matrix, VB8[2]))

## 48 buckets per day ≡ 10 min in Calendar time
VB48 = MakeVolumeTimeData(prices, volume, 48)

MMVB48 = NUFFTcorrDKFGG(convert(Matrix, VB48[1]), convert(Matrix, VB48[2]))[1]
HYVB48 = HYcorrFull(convert(Matrix, VB48[1]), convert(Matrix, VB48[2]))

## 480 buckets per day ≡ 1 min in Calendar time
VB480 = MakeVolumeTimeData(prices, volume, 480)

MMVB480 = NUFFTcorrDKFGG(convert(Matrix, VB480[1]), convert(Matrix, VB480[2]))[1]
HYVB480 = HYcorrFull(convert(Matrix, VB480[1]), convert(Matrix, VB480[2]))

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

## Volume Time:
# 8 buckets
plot(MMVB8, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :orange, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho^{ij}|} = $(round(summarise(MMVB8)[1], digits = 4)) \\pm $(round(summarise(MMVB8)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/VB8.png")

# 48 buckets
plot(MMVB48, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :orange, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho^{ij}|} = $(round(summarise(MMVB48)[1], digits = 4)) \\pm $(round(summarise(MMVB48)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/VB48.png")

# 480 buckets
plot(MMVB480, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :orange, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho^{ij}|} = $(round(summarise(MMVB480)[1], digits = 4)) \\pm $(round(summarise(MMVB480)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/VB480.png")
