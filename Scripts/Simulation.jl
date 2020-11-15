## Author: Patrick Chang
# Script file to investigate the Epps effect in event time
# using a Hawkes process to simulate the prices.
# We use the price model by Bacry et al. (2013) and compare
# the MM and HY estimator.

## Preamble

using Pipe, StatsBase, ProgressMeter, Distributions, JLD

cd("/Users/patrickchang1/PCRBTG-VT")

# Estimators
include("../Functions/Correlation Estimators/Dirichlet/NUFFTcorrDK-FGG.jl")
include("../Functions/Correlation Estimators/HY/HYcorr.jl")

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

# Helpful functions
#------------------------
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

# Making the prices
#------------------------
# Construct asynchronously sampled prices for MM and HY in calendar time
function CalendarTimeMMHYPrices(Asset1, Asset2, P01, P02)
    # Get trading times and prices from each asset
    t1 = Asset1[2]; t2 = Asset2[2]
    P1 = Asset1[1]; P2 = Asset2[1]
    # Combine all trading times
    times = [0;t1;t2]
    sort!(times)
    # Initialize variables
    n = length(times)
    P = fill(NaN, n, 2)
    tt = fill(NaN, n, 2)
    # Starting values
    P[1, 1] = P01; tt[1, 1] = 0
    P[1, 2] = P02; tt[1, 2] = 0
    # Counters
    i_1 = 1
    i_2 = 1
    # Loop through each event
    for i in 2:n
        if in(times[i], t1)
            P[i,1] = P1[i_1]
            tt[i,1] = t1[i_1]
            i_1 += 1
        else
            P[i,2]= P2[i_2]
            tt[i,2]= t2[i_2]
            i_2 += 1
        end
    end
    return P, tt
end
# Construct uniformly sampled prices for RV in calendar time
# apply PTI to deal with asynchrony
function CalendarTimeRVPrices(P0, τ, T, uptick, downtick)
    t = collect(0:τ:T)
    n = length(t)
    P = zeros(n,1)

    for i in 1:n
        up = findall(uptick .<= t[i])
        down = findall(downtick .<= t[i])
        P[i] = P0 + length(up) - length(down)
    end
    return P
end
# Construct asynchronously sampled prices for MM and HY in event time
function EventTimeMMHYPrices(Asset1, Asset2, P01, P02)
    # Get trading times and prices from each asset
    t1 = Asset1[2]; t2 = Asset2[2]
    P1 = Asset1[1]; P2 = Asset2[1]
    # Combine all trading times
    times = [0;t1;t2]
    sort!(times)
    # Initialize variables
    n = length(times)
    P = fill(NaN, n, 2)
    tt = fill(NaN, n, 2)
    # Starting values
    P[1, 1] = P01; tt[1, 1] = 1
    P[1, 2] = P02; tt[1, 2] = 1
    # Counters
    i_1 = 1
    i_2 = 1
    # Loop through each event
    for i in 2:n
        if in(times[i], t1)
            P[i,1] = P1[i_1]
            tt[i,1] = i
            i_1 += 1
        else
            P[i,2]= P2[i_2]
            tt[i,2]= i
            i_2 += 1
        end
    end
    return P, tt
end
# Construct uniformly sampled prices for RV in  event time
# apply PTI to deal with asynchrony
function EventTimeRVPrice(Asset1, Asset2, P01, P02)
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
function VolumeTimePrice(Asset1, Asset2, num::Int, seed)
    ### num: the number of samples we want
    # Get prices from each asset
    P1 = exp.(Asset1[1]); P2 = exp.(Asset2[1])
    # Simulate random volumes to a max of maxvol
    Random.seed!(seed)
    V1 = round.(RPareto(20, 1.7, length(P1)))
    Random.seed!(seed + length(P1))
    V2 = round.(RPareto(20, 1.7, length(P2)))
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

# Function to extract estimates under calendar time over different sampling intervals
function getCalendarTimeCorrelations(Hawkesdata, Asyndata, dt, T)
    n = length(dt)
    # Initialize storage variables for the estimators
    MM = zeros(n, 1)
    RVest = zeros(n, 1)

    for i in 1:length(dt)
        # MM time scale
        N = Int(floor((T/dt[i] - 1)/2))
        # RV data
        p1 = CalendarTimeRVPrices(0, dt[i], T, Hawkesdata[1], Hawkesdata[2])
        p2 = CalendarTimeRVPrices(0, dt[i], T, Hawkesdata[3], Hawkesdata[4])
        tempdata = exp.([p1 p2])
        # tempdata = kskip(Syndata, dt[i])
        # Compute the results
        MM[i] = NUFFTcorrDKFGG(exp.(Asyndata[1]), Asyndata[2], N = N)[1][1,2]
        RVest[i] = RV(tempdata)[1][1,2]
    end
    HY = HYcorr(exp.(Asyndata[1][:,1]), exp.(Asyndata[1][:,2]), Asyndata[2][:,1], Asyndata[2][:,2])[1][1,2]
    return MM, RVest, HY
end

# Function to extract estimates under event time over different sampling intervals
function getEventTimeCorrelations(Syndata, Asyndata, dt)
    # Get range of time scales
    prices = Syndata[1]
    t = Syndata[2]
    n = length(dt)
    # Initialize storage variables for the estimators
    MM = zeros(n, 1)
    RVest = zeros(n, 1)
    # Loop through the ranges
    for i in 1:n
        # MM time scale
        N = Int(floor((size(Asyndata[1])[1]/dt[i] - 1)/2))
        # HY k-skip
        tempdata = kskip(prices, dt[i])
        # Compute the results
        MM[i] = NUFFTcorrDKFGG(exp.(Asyndata[1]), Asyndata[2], N = N)[1][1,2]
        RVest[i] = RV(exp.(tempdata))[1][1,2]
    end
    HY = HYcorr(exp.(Asyndata[1][:,1]), exp.(Asyndata[1][:,2]), Asyndata[2][:,1], Asyndata[2][:,2])[1][1,2]
    return MM, RVest, HY
end

# Function to extract estimates under volume time over different bucket sizes
function getVolumeTimeCorrelations(Asset1, Asset2, dt, T, seed)
    # Get range of time scales
    n = length(dt)
    # Initialize storage variables for the estimators
    MM = zeros(n, 1)
    RVest = zeros(n, 1)
    HY = zeros(n, 1)
    # Loop through the ranges
    for i in 1:n
        # Create volume data for given certain amount of samples
        tempdata = VolumeTimePrice(Asset1, Asset2, Int(floor(T/dt[i])), seed)
        t = [collect(1:1:size(tempdata)[1]) collect(1:1:size(tempdata)[1])]
        # Compute the results
        # MM[i] = NUFFTcorrDKFGG(exp.(tempdata), t)[1][1,2]
        RVest[i] = RV(exp.(tempdata))[1][1,2]
        # HY[i] = HYcorr(exp.(tempdata[:,1]), exp.(tempdata[:,2]), t[:,1], t[:,2])[1][1,2]
    end
    # return MM, RVest, HY
    return RVest
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
    MMET = zeros(reps, length(domain)); MMVT = zeros(reps, length(domain)); MMCT = zeros(reps, length(domain))
    RVET = zeros(reps, length(domain)); RVVT = zeros(reps, length(domain)); RVCT = zeros(reps, length(domain))
    HYET = zeros(reps, 1); HYCT = zeros(reps, 1); HYVT = zeros(reps, length(domain))
    # Loop through the different cases
    @showprogress "Computing..." for i in 1:reps
        # Simulate the processs
        t = simulateHawkes(lambda0_1, alpha_1, beta_1, T, seed = seeds[i])
        # Create asynchronous calendar time data
        CTMMHYData = CalendarTimeMMHYPrices(makePrices(0, t[1], t[2]), makePrices(0, t[3], t[4]), 0, 0)
        # Create synchronous event time data
        ETRVData = EventTimeRVPrice(makePrices(0, t[1], t[2]), makePrices(0, t[3], t[4]), 0, 0)
        # Creat asynchronous event time data
        ETMMHYData = EventTimeMMHYPrices(makePrices(0, t[1], t[2]), makePrices(0, t[3], t[4]), 0, 0)
        # Obtain the results
        CalendarTimeRes = getCalendarTimeCorrelations(t, CTMMHYData, domain, T)
        EventTimeRes = getEventTimeCorrelations(ETRVData, ETMMHYData, domain)
        VolumeTimeRes = getVolumeTimeCorrelations(makePrices(0, t[1], t[2]), makePrices(0, t[3], t[4]), domain, T)
        # Store the results
        MMCT[i,:] = CalendarTimeRes[1]; RVCT[i,:] = CalendarTimeRes[2]; HYCT[i] = CalendarTimeRes[3]
        MMET[i,:] = EventTimeRes[1]; RVET[i,:] = EventTimeRes[2]; HYET[i] = EventTimeRes[3]
        MMVT[i,:] = VolumeTimeRes[1]; RVVT[i,:] = VolumeTimeRes[2]; HYVT[i,:] = VolumeTimeRes[3]
    end
    return MMCT, RVCT, HYCT, MMET, RVET, HYET, MMVT, RVVT, HYVT
end

function getRes(reps)
    # Set up parameters
    T = 3600*20
    par_1 = BarcyParams(0.015, 0.023, 0.05, 0.11)
    lambda0_1 = par_1[1]; alpha_1 = par_1[2]; beta_1 = par_1[3]
    domain = collect(1:1:100)
    Random.seed!(2020)
    seeds = Int.(floor.(rand(reps) .* 1000000))
    # Initialize storage variables
    RVVT = zeros(reps, length(domain))
    # Loop through the different cases
    @showprogress "Computing..." for i in 1:reps
        # Simulate the processs
        t = simulateHawkes(lambda0_1, alpha_1, beta_1, T, seed = seeds[i])
        # Obtain the results
        VolumeTimeRes = getVolumeTimeCorrelations(makePrices(0, t[1], t[2]), makePrices(0, t[3], t[4]), domain, T, seeds[i])
        # Store the results
        RVVT[i,:] = VolumeTimeRes
    end
    return RVVT
end
res2 = getRes(reps)

test = []
for i in 1:6
    push!(test, res[i])
end
push!(test, res2)
push!(test, res2)
push!(test, res2)

res = test

## Obtain the results
reps = 100
res = getRes(reps)

# Save results
save("Computed Data/Simulation.jld", "res", res)
# Load results
res = load("Computed Data/Simulation.jld")
res = res["res"]

# Theoretical Epps effect from Hawkes process in Calendar time
dt = collect(1:1:100)
theoretical = zeros(length(dt), 1)
for i in 1:length(dt)
    theoretical[i] = theoreticalEpps(dt[i], 0.015, 0.023, 0.05, 0.11)
end

## Plot the results
ρ = theoreticalCorr(0.023, 0.05, 0.11)
q = quantile.(TDist(reps-1), [0.975])

# Calendar time
p1 = plot(dt, mean(res[1], dims=1)', ribbon = (q .* std(res[1], dims=1)'), fillalpha = .3, legend = :bottomright, color = :blue, line=(1, [:dash]), label = L"\textrm{MM}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)), dpi = 300)
plot!(p1, dt, mean(res[2], dims=1)', ribbon = (q .* std(res[2], dims=1)'), fillalpha=.3, color = :red, line=(1, [:dash]), label = L"\textrm{RV}", marker = ([:x :d], 3, 0.8, stroke(3, :red)))
plot!(p1, dt, theoretical, color = :black, line=(2, [:solid]), label = L"\textrm{Theoretical Epps effect}")
hline!(p1, [mean(res[3])], ribbon = (q .* std(res[3], dims=1)'), fillalpha=.3, color = :brown, line=(1, [:dash]), label = L"\textrm{HY}")
hline!(p1, [ρ], color = :black, line=(2, [:dot]), label = L"\textrm{Limiting } \rho")
xlabel!(p1, L"\textrm{Sampling interval (calendar time)}")
ylabel!(p1, L"\tilde{\rho}_{\Delta t}^{ij}")
savefig(p1, "Plots/SimCT.svg")

# Transaction time
p2 = plot(dt, mean(res[4], dims=1)', ribbon = (q .* std(res[4], dims=1)'), fillalpha = .3, legend = :bottomright, color = :blue, line=(1, [:dash]), label = L"\textrm{MM}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)), dpi = 300)
plot!(p2, dt, mean(res[5], dims=1)', ribbon = (q .* std(res[5], dims=1)'), fillalpha=.3, color = :red, line=(1, [:dash]), label = L"\textrm{RV}", marker = ([:x :d], 3, 0.8, stroke(3, :red)))
hline!(p2, [mean(res[6])], ribbon = (q .* std(res[6], dims=1)'), fillalpha=.3, color = :brown, line=(1, [:dash]), label = L"\textrm{HY}")
hline!(p2, [ρ], color = :black, line=(2, [:dot]), label = L"\textrm{Limiting } \rho")
xlabel!(p2, L"\textrm{Sampling interval (event time)}")
ylabel!(p2, L"\tilde{\rho}_{\Delta t}^{ij}")
savefig(p2, "Plots/SimET.svg")

# Volume time
p3 = plot(dt, mean(res[7], dims=1)', ribbon = (q .* std(res[7], dims=1)'), fillalpha = .3, legend = :topright, color = :blue, line=(1, [:dash]), label = L"\textrm{MM}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)), dpi = 300)
plot!(p3, dt, mean(res[8], dims=1)', ribbon = (q .* std(res[8], dims=1)'), fillalpha=.3, color = :red, line=(1, [:dash]), label = L"\textrm{RV}", marker = ([:x :d], 3, 0.8, stroke(3, :red)))
plot!(p3, dt, mean(res[9], dims=1)', ribbon = (q .* std(res[9], dims=1)'), fillalpha=.3, color = :brown, line=(1, [:dash]), label = L"\textrm{HY}", marker = ([:x :d], 3, 0.8, stroke(3, :red)))
hline!(p3, [ρ], color = :black, line=(2, [:dot]), label = L"\textrm{Limiting } \rho")
xlabel!(p3, L"\textrm{Sampling interval (volume time)}")
ylabel!(p3, L"\tilde{\rho}_{\Delta t}^{ij}")
savefig(p3, "Plots/SimVT.svg")

# Comparison
p4 = plot(dt, mean(res[1], dims=1)', legend = :right, color = :blue, line=(1, [:dash]), label = L"\textrm{CT MM}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)), dpi = 300)
plot!(p4, dt, mean(res[2], dims=1)', color = :red, line=(1, [:dash]), label = L"\textrm{CT RV}", marker = ([:x :d], 3, 0.8, stroke(3, :red)))
plot!(p4, dt, mean(res[4], dims=1)', color = :green, line=(1, [:dash]), label = L"\textrm{ET MM}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)), dpi = 300)
plot!(p4, dt, mean(res[5], dims=1)', color = :orange, line=(1, [:dash]), label = L"\textrm{ET RV}", marker = ([:x :d], 3, 0.8, stroke(3, :red)))
plot!(p4, dt, mean(res[8], dims=1)', color = :pink, line=(1, [:dash]), label = L"\textrm{VT RV}", marker = ([:+ :d], 3, 0.8, stroke(3, :red)))
hline!(p4, [mean(res[3])], color = :brown, line=(1, [:dash]), label = L"\textrm{CT HY}")
hline!(p4, [mean(res[6])], color = :lightblue, line=(2, [:dashdot]), label = L"\textrm{ET HY}")
plot!(p4, dt, theoretical, color = :black, line=(2, [:solid]), label = L"\textrm{Theoretical Epps effect}")
hline!(p4, [ρ], color = :black, line=(2, [:dot]), label = L"\textrm{Limiting } \rho")
xlabel!(p4, L"\textrm{Sampling interval}")
ylabel!(p4, L"\tilde{\rho}_{\Delta t}^{ij}")
savefig(p4, "Plots/SimComp.svg")


## Ploting the price paths

par = BarcyParams(0.015, 0.023, 0.05, 0.11)
lambda0 = par[1]; alpha = par[2]; beta = par[3]
T = 300
t = simulateHawkes(lambda0, alpha, beta, T)

## Calendar time
# PIT
CTRVp1 = CalendarTimeRVPrices(0, 10, T, t[1], t[2])
CTRVp2 = CalendarTimeRVPrices(0, 10, T, t[3], t[4])
CTRVPrices = [CTRVp1 CTRVp2]
CTRVtimes = [collect(0:10:T) collect(0:10:T)]

plot(CTRVtimes, CTRVPrices, seriestype = :steppost, color = [:blue :red], marker = :dot, legend = :topleft, label = [L"X^1" L"X^2"])
xlabel!(L"\textrm{Calendar time}")
ylabel!(L"X^i_t")
savefig("Plots/CalendarTimeRVPrices.svg")

# Asyn
CTMMHY = CalendarTimeMMHYPrices(makePrices(0, t[1], t[2]), makePrices(0, t[3], t[4]), 0, 0)

plot(filter(!isnan, CTMMHY[2][:,1]), filter(!isnan, CTMMHY[1][:,1]), seriestype = :steppost, color = :blue, marker = :dot, legend = :topleft, label = L"X^1")
plot!(filter(!isnan, CTMMHY[2][:,2]), filter(!isnan, CTMMHY[1][:,2]), seriestype = :steppost, color = :red, marker = :dot, label = L"X^2")
xlabel!(L"\textrm{Calendar time}")
ylabel!(L"X^i_t")
savefig("Plots/CalendarTimeMMHYPrices.svg")

## Event time
# PIT
ETRV = EventTimeRVPrice(makePrices(0, t[1], t[2]), makePrices(0, t[3], t[4]), 0, 0)

plot(ETRV[2], ETRV[1], seriestype = :steppost, linetype = :steppost, marker = :dot, color=[:blue :red], legend = :topleft, label = [L"X^1" L"X^2"])
xlabel!(L"\textrm{Event time}")
ylabel!(L"X^i_t")
savefig("Plots/EventTimeRVPrices.svg")

# Asyn
ETMMHY = EventTimeMMHYPrices(makePrices(0, t[1], t[2]), makePrices(0, t[3], t[4]), 0, 0)

plot(filter(!isnan, ETMMHY[2][:,1]), filter(!isnan, ETMMHY[1][:,1]), seriestype = :steppost, color = :blue, marker = :dot, legend = :topleft, label = L"X^1")
plot!(filter(!isnan, ETMMHY[2][:,2]), filter(!isnan, ETMMHY[1][:,2]), seriestype = :steppost, color = :red, marker = :dot, label = L"X^2")
xlabel!(L"\textrm{Event time}")
ylabel!(L"X^i_t")
savefig("Plots/EventTimeMMHYPrices.svg")

## Volume time

VT30 = VolumeTimePrice(makePrices(0, t[1], t[2]), makePrices(0, t[3], t[4]), 30, 1)

plot(1:30, VT30, seriestype = :steppost, linetype = :steppost, marker = :dot, color=[:blue :red], legend = :topleft, label = [L"X^1" L"X^2"], xlims = (0,31))
xlabel!(L"\textrm{Volume time}")
ylabel!(L"X^i_t")
savefig("Plots/VolumeTime30Prices.svg")


VT60 = VolumeTimePrice(makePrices(0, t[1], t[2]), makePrices(0, t[3], t[4]), 60, 1)

plot(1:60, VT60, seriestype = :steppost, linetype = :steppost, marker = :dot, color=[:blue :red], legend = :topleft, label = [L"X^1" L"X^2"])
xlabel!(L"\textrm{Volume time}")
ylabel!(L"X^i_t")
savefig("Plots/VolumeTime60Prices.svg")
