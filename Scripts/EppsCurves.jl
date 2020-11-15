## Author: Patrick Chang
# Script file to investigate the Epps curves for different
# methods of aggregating data. Consider MM calendar time,
# PTI in calendar time, event time and volume time. Also consider HY as baseline.

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
# Combine the data
data = Dict()
push!(data, :prices => prices); push!(data, :times => times); push!(data, :volume => volume)

## Useful functions

# Function to create the volume time data of with a specific number of buckets per day
function MakeVolumeTimeData(prices::DataFrame, volume::DataFrame, numbuckets::Integer)
    # Create the new dataframe to store the bar data
    prices_VT = zeros(numbuckets*5, size(prices)[2])
    # Loop through each column
    for i in 1:size(prices)[2]
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
    for i in 2:n
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
function kskipEvents(data, numevents::Int64)
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

# Function to create the 1 second interval PTI prices
function MakePTI(prices::DataFrame, times::DataFrame)
    # Create the new dataframe to store the bar data
    master_df_prices = DataFrame()
    # Loop through each column
    for i in 1:size(prices)[2]
        # Extract temporary data
        temp_price = filter!(!isnan, prices[:,i])
        temp_time = filter!(!isnan, times[:,i])
        # Convert unix time to date time
        dates = Date.(unix2datetime.(temp_time))
        dates_unique = unique(dates)
        # Initialise temporary vectors
        temp_df_prices = []
        temp_df_times = []
        # Loop through each day
        for j in 1:length(dates_unique)
            # Create extract data from each day
            tempday = dates_unique[j]
            inds = findall(x -> x == tempday, dates)
            tempdata_price = temp_price[inds]
            tempdata_time = unix2datetime.(temp_time)[inds]
            # Only keep data within continuous trading
            start = DateTime(tempday) + Hour(9-2)
            close = DateTime(tempday) + Hour(16-2) + Minute(50)
            # Create the bars
            bars = collect(start:Second(1):close)
            # Cheap hack. ASSUME: that events all start together, achieved by bringing the first obs in each ticker to the start
            tempdata_price = [tempdata_price[1]; tempdata_price]
            tempdata_time = [start; tempdata_time]
            # Loop through each bar
            for k in 2:length(bars)
                push!(temp_df_prices, tempdata_price[findall(x-> x<=bars[k], tempdata_time)[end]])
            end
        end
        # Store the data into master_df
        master_df_prices[:,i] = temp_df_prices
    end
    return master_df_prices
end

# k-skip sampling
function kskipPTI(data, Δt)
    # Get dimensions
    (n, m) = size(data)
    # Get the indicies
    inds = collect(1:Δt:n)
    # Extract the prices
    P = data[inds,:]
    return P
end

## Function to streamline everything together
function getCurves(A1::Symbol, A2::Symbol, data::Dict)
    # Extract the tickers from the data
    P = data[:prices][:,[A1; A2]]; t = data[:times][:,[A1; A2]]; V = data[:volume][:,[A1; A2]]
    # Set the range for investigation
    dt = collect(5:5:600)
    # Initialize the storage variables
    MMcalendar = zeros(length(dt), 1); RVcalendar = zeros(length(dt), 1)
    MMET = zeros(length(dt), 1); MMVT = zeros(length(dt), 1)
    RVET = zeros(length(dt), 1); HY = zeros(1, 1)
    # Create event time data and calendar time PTI
    EventTime = MakeEventTimeData(P)
    CalendarTime = MakePTI(P, t)
    # Loop through the various time scales
    @showprogress "Computing..." for i in 1:length(dt)
        # Get time scales for MM estimator in calendar and event time
        Ncal = Int(floor((141000/dt[i] - 1)/2))
        Nevent = Int(floor((size(EventTime[1])[1]/dt[i] - 1)/2))
        # Get calendar and volume time data
        CTk = kskipPTI(CalendarTime, dt[i])
        VT = MakeVolumeTimeData(P, V, Int(floor(((60*7+50)*60)/dt[i])))
        ETk = kskipPTI(EventTime[1], dt[i])
        # Extract results
        MMcalendar[i] = NUFFTcorrDKFGG(convert(Matrix, P), convert(Matrix, t), N = Ncal)[1][1,2]
        RVcalendar[i] = RV(convert(Matrix, CTk))[1][1,2]
        MMET[i] = NUFFTcorrDKFGG(convert(Matrix, EventTime[1]), convert(Matrix, EventTime[2]), N = Nevent)[1][1,2]
        RVET[i] = RV(convert(Matrix, ETk))[1][1,2]
        MMVT[i] = NUFFTcorrDKFGG(VT[1], VT[2])[1][1,2]
    end
    HY = HYcorrFull(convert(Matrix, P), convert(Matrix, t))[1,2]
    return MMcalendar, RVcalendar, MMET, RVET, MMVT, HY
end



dt = collect(5:5:600)

SBKFSR = getCurves(:FSR, :SBK, data)

p1 = plot(dt, SBKFSR[1], legend = :outertopright, color = :blue, line=(1, [:solid]), label = L"\textrm{MM Calendar}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)), dpi = 300, size = (800, 500))
plot!(p1, dt, SBKFSR[2], color = :red, line=(1, [:solid]), label = L"\textrm{RV Calendar}", marker = ([:+ :d], 3, 0.8, stroke(3, :red)))
plot!(p1, dt, SBKFSR[3], color = :blue, line=(1, [:dash]), label = L"\textrm{MM Event time}", marker = ([:x :d], 3, 0.8, stroke(3, :blue)))
plot!(p1, dt, SBKFSR[4], color = :red, line=(1, [:dash]), label = L"\textrm{RV Event time}", marker = ([:x :d], 3, 0.8, stroke(3, :red)))
plot!(p1, dt, SBKFSR[5], color = :green, line=(1, [:dot]), label = L"\textrm{Volume time}", marker = ([:circ :d], 3, 0.8, stroke(3, :green)))
hline!(p1, [SBKFSR[6]], color = :brown, line=(1, [:dash]), label = L"\textrm{HY}")


NEGABG = getCurves(:NED, :ABG, data)

p2 = plot(dt, NEGABG[1], legend = :outertopright, color = :blue, line=(1, [:solid]), label = L"\textrm{MM Calendar}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)), dpi = 300, size = (800, 500))
plot!(p2, dt, NEGABG[2], color = :red, line=(1, [:solid]), label = L"\textrm{RV Calendar}", marker = ([:+ :d], 3, 0.8, stroke(3, :red)))
plot!(p2, dt, NEGABG[3], color = :blue, line=(1, [:dash]), label = L"\textrm{MM Event time}", marker = ([:x :d], 3, 0.8, stroke(3, :blue)))
plot!(p2, dt, NEGABG[4], color = :red, line=(1, [:dash]), label = L"\textrm{RV Event time}", marker = ([:x :d], 3, 0.8, stroke(3, :red)))
plot!(p2, dt, NEGABG[5], color = :green, line=(1, [:dot]), label = L"\textrm{Volume time}", marker = ([:circ :d], 3, 0.8, stroke(3, :green)))
hline!(p2, [NEGABG[6]], color = :brown, line=(1, [:dash]), label = L"\textrm{HY}")
