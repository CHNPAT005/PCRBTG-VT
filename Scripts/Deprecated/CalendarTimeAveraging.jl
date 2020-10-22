## Author: Patrick Chang
# Script file to investigate the correlation dynamic for 10 JSE equity
# using calendar time averaging. We consider the closing and VWAP price
# for 1Hr, 10min, 1min and TAQ data.

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

## Function to make bar data
#---------------------------------------------------------------------------

# Function to create the OHLCV + VWAP bar data using transaction data
# Note that bar size only takes in minutes as the argument
function MakeTransactionBars(prices::DataFrame, times::DataFrame, volume::DataFrame, barsize::Integer)
    # Create the new dataframe to store the bar data
    master_df_prices_closing = DataFrame()
    master_df_prices_VWAP = DataFrame()
    master_df_times = DataFrame()
    # Loop through each column
    @showprogress "Computing..." for i in 1:size(prices)[2]
        # Extract temporary data
        temp_price = filter!(!isnan, prices[:,i])
        temp_time = filter!(!isnan, times[:,i])
        temp_vol = filter!(!isnan, volume[:,i])
        # Convert unix time to date time
        dates = Date.(unix2datetime.(temp_time))
        dates_unique = unique(dates)
        # Initialise temporary vectors
        temp_df_prices_closing = []
        temp_df_prices_VWAP = []
        temp_df_times = []
        # Loop through each day
        for j in 1:length(dates_unique)
            # Create extract data from each day
            tempday = dates_unique[j]
            inds = findall(x -> x == tempday, dates)
            tempdata_price = temp_price[inds]
            tempdata_time = unix2datetime.(temp_time)[inds]
            tempdata_vol = temp_vol[inds]
            # Only keep data within continuous trading
            start = DateTime(tempday) + Hour(9-2)
            close = DateTime(tempday) + Hour(16-2) + Minute(50)
            # Create the bars
            bars = collect(start:Minute(barsize):close)
            # Loop through each bar
            for k in 2:length(bars)
                # filter out data within the bar
                bardata_price = tempdata_price[findall(x-> x>bars[k-1] && x<=bars[k], tempdata_time)]
                bardata_vol = tempdata_vol[findall(x-> x>bars[k-1] && x<=bars[k], tempdata_time)]
                # Check if bar is empty
                if !isempty(bardata_price)
                    # Not empty; make Closing
                    closing = bardata_price[end]
                    temp_df_prices_closing = [temp_df_prices_closing; closing]
                    # make VWAP
                    VWAP = sum(bardata_price .* bardata_vol) / sum(bardata_vol)
                    temp_df_prices_VWAP = [temp_df_prices_VWAP; VWAP]
                    # add the times
                    temp_df_times = [temp_df_times; datetime2unix(bars[k])]
                else
                    # Add NaN
                    temp_df_prices_closing = [temp_df_prices_closing; NaN]
                    temp_df_prices_VWAP = [temp_df_prices_VWAP; NaN]
                    temp_df_times = [temp_df_times; NaN]
                end
            end
        end
        # Store the data into master_df
        master_df_prices_closing[:,i] = temp_df_prices_closing
        master_df_prices_VWAP[:,i] = temp_df_prices_VWAP
        master_df_times[:,i] = temp_df_times
    end
    return master_df_prices_closing, master_df_prices_VWAP, master_df_times
end

## 1 Hour bar:
bar1hr = MakeTransactionBars(prices, times, volume, 60)

# Closing
MMClosing1hr = NUFFTcorrDKFGG(convert(Matrix,bar1hr[1]), convert(Matrix,bar1hr[3]))[1]
HYClosing1hr = HYcorrFull(convert(Matrix,bar1hr[1]), convert(Matrix,bar1hr[3]))

# VWAP
MMVWAP1hr = NUFFTcorrDKFGG(convert(Matrix,bar1hr[2]), convert(Matrix,bar1hr[3]))[1]
HYVWAP1hr = HYcorrFull(convert(Matrix,bar1hr[2]), convert(Matrix,bar1hr[3]))


## 10 min bar:
bar10min = MakeTransactionBars(prices, times, volume, 10)

# Closing
MMClosing10min = NUFFTcorrDKFGG(convert(Matrix,bar10min[1]), convert(Matrix,bar10min[3]))[1]
HYClosing10min = HYcorrFull(convert(Matrix,bar10min[1]), convert(Matrix,bar10min[3]))

# VWAP
MMVWAP10min = NUFFTcorrDKFGG(convert(Matrix,bar10min[2]), convert(Matrix,bar10min[3]))[1]
HYVWAP10min = HYcorrFull(convert(Matrix,bar10min[2]), convert(Matrix,bar10min[3]))


## 1 Hour bar:
bar1min = MakeTransactionBars(prices, times, volume, 1)

# Closing
MMClosing1min = NUFFTcorrDKFGG(convert(Matrix,bar1min[1]), convert(Matrix,bar1min[3]))[1]
HYClosing1min = HYcorrFull(convert(Matrix,bar1min[1]), convert(Matrix,bar1min[3]))

# VWAP
MMVWAP1min = NUFFTcorrDKFGG(convert(Matrix,bar1min[2]), convert(Matrix,bar1min[3]))[1]
HYVWAP1min = HYcorrFull(convert(Matrix,bar1min[2]), convert(Matrix,bar1min[3]))

## TAQ:

# TAQ
MMTAQ = NUFFTcorrDKFGG(convert(Matrix,prices), convert(Matrix,times))[1]
HYTAQ = HYcorrFull(convert(Matrix,prices), convert(Matrix,times))

# ## Save the results
# save("Computed Data/CalendarTime.jld",
# "MMClosing1hr", MMClosing1hr, "HYClosing1hr", HYClosing1hr, "MMVWAP1hr", MMVWAP1hr, "HYVWAP1hr", HYVWAP1hr,
# "MMClosing10min", MMClosing10min. "HYClosing10min", HYClosing10min, "MMVWAP10min", MMVWAP10min, "HYVWAP10min", HYVWAP10min,
# "MMClosing1min", MMClosing1min, "HYClosing1min", HYClosing1min, "MMVWAP1min", MMVWAP1min, "HYVWAP1min", HYVWAP1min,
# "MMTAQ", MMTAQ, "HYAQ", HYAQ)
#
# ## Load the results
# res = load("Computed Data/CalendarTime.jld")
# MMClosing1hr = res["MMClosing1hr"]; HYClosing1hr = res["HYClosing1hr"]
# MMVWAP1hr = res["MMVWAP1hr"]; HYVWAP1hr = res["HYVWAP1hr"]
# MMClosing10min = res["MMClosing10min"]; HYClosing10min = res["HYClosing10min"]
# MMVWAP10min = res["MMVWAP10min"]; HYVWAP10min = res["HYVWAP10min"]
# MMClosing1min = res["MMClosing1min"]; HYClosing1min = res["HYClosing1min"]
# MMVWAP1min = res["MMVWAP1min"]; HYVWAP1min = res["HYVWAP1min"]
# MMTAQ = res["MMTAQ"]; HYAQ = res["HYAQ"]


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

## Closing:
# MM
plot(MMClosing1hr, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho_{ij}|} = $(round(summarise(MMClosing1hr)[1], digits = 4)) \\pm $(round(summarise(MMClosing1hr)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/MMClosing1hr.png")

plot(MMClosing10min, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho_{ij}|} = $(round(summarise(MMClosing10min)[1], digits = 4)) \\pm $(round(summarise(MMClosing10min)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/MMClosing10min.png")

plot(MMClosing1min, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho_{ij}|} = $(round(summarise(MMClosing1min)[1], digits = 4)) \\pm $(round(summarise(MMClosing1min)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/MMClosing1min.png")

# HY
plot(HYClosing1hr, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho_{ij}|} = $(round(summarise(HYClosing1hr)[1], digits = 4)) \\pm $(round(summarise(HYClosing1hr)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/HYClosing1hr.png")

plot(HYClosing10min, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho_{ij}|} = $(round(summarise(HYClosing10min)[1], digits = 4)) \\pm $(round(summarise(HYClosing10min)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/HYClosing10min.png")

plot(HYClosing1min, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho_{ij}|} = $(round(summarise(HYClosing1min)[1], digits = 4)) \\pm $(round(summarise(HYClosing1min)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/HYClosing1min.png")

## VWAP
# MM
plot(MMVWAP1hr, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho_{ij}|} = $(round(summarise(MMVWAP1hr)[1], digits = 4)) \\pm $(round(summarise(MMVWAP1hr)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/MMVWAP1hr.png")

plot(MMVWAP10min, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho_{ij}|} = $(round(summarise(MMVWAP10min)[1], digits = 4)) \\pm $(round(summarise(MMVWAP10min)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/MMVWAP10min.png")

plot(MMVWAP1min, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho_{ij}|} = $(round(summarise(MMVWAP1min)[1], digits = 4)) \\pm $(round(summarise(MMVWAP1min)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/MMVWAP1min.png")

# HY
plot(HYVWAP1hr, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho_{ij}|} = $(round(summarise(HYVWAP1hr)[1], digits = 4)) \\pm $(round(summarise(HYVWAP1hr)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/HYVWAP1hr.png")

plot(HYVWAP10min, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho_{ij}|} = $(round(summarise(HYVWAP10min)[1], digits = 4)) \\pm $(round(summarise(HYVWAP10min)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/HYVWAP10min.png")

plot(HYVWAP1min, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho_{ij}|} = $(round(summarise(HYVWAP1min)[1], digits = 4)) \\pm $(round(summarise(HYVWAP1min)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/HYVWAP1min.png")

## TAQ
# MM
plot(MMTAQ, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho_{ij}|} = $(round(summarise(MMTAQ)[1], digits = 4)) \\pm $(round(summarise(MMTAQ)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/MMTAQ.png")

# HY
plot(HYTAQ, st=:heatmap, clim=(-1,1), color=cgrad([:blue, :green, :yellow, :red]), colorbar_title=L"\rho^{ij}", xticks = (1:10, tickers), yticks = (1:10, tickers), dpi = 300, size = (800, 700), tickfontsize = 15)
plot!(annotations=(5, 4, Plots.text(latexstring("\$\\overline{|\\rho_{ij}|} = $(round(summarise(HYAQ)[1], digits = 4)) \\pm $(round(summarise(HYAQ)[2], digits = 4)) \$"), :left, 20)))
# savefig("Plots/HYTAQ.png")
