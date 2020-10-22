## Author: Patrick Chang
# Script file to investigate the effect of
# 1) Volatility clustering using a GARCH(1,1), and
# 2) Mean reversion using a Ornstein Uhlenbeck process
# Here we look at both the synchronous and asynchronous case
# with the asynchronous case induced by down-sampling 20%.

## Preamble

using Plots, Random, ProgressMeter, StatsBase, JLD, LaTeXStrings

cd("/Users/patrickchang1/PCRBTG-VT")

# Estimators
include("../Functions/Correlation Estimators/Dirichlet/NUFFTcorrDK-FGG.jl")
include("../Functions/Correlation Estimators/HY/HYcorr.jl")

# SDE
include("../Functions/SDEs/Garch.jl")
include("../Functions/SDEs/OU.jl")

## Fucntion to strealine the process
# The function uses a GARCH(1,1) process.
# Asynchrony is achieved by down-sampling by 20%.
function VolClust(num)
    # Setup
    cor = range(-0.99, 0.99, length = num)
    MM = zeros(num, 1)
    HY = zeros(num, 1)
    MMasyn = zeros(num, 1)
    HYasyn = zeros(num, 1)
    # Parameters
    θ = [0.035, 0.054]
    λ = [0.296, 0.48]
    ω = [0.636, 0.476]
    seed = 1:num
    # Run through each case
    @showprogress "Computing..." for i in 1:num
        # Simulate the GARCH
        p = Garch(10000, θ, λ, ω, cor[i]; seed = seed[i])
        t = [collect(1:1:10000.0) collect(1:1:10000.0)]
        # Compute the correlation
        MM[i] = NUFFTcorrDKFGG(p, t)[1][1,2]
        HY[i] = HYcorr(p[:,1], p[:,2], t[:,1], t[:,2])[1][1,2]
        # Downsample
        Random.seed!(seed[i])
        rm1 = sample(1:10000, Int(10000*0.2))
        Random.seed!(seed[i] + num)
        rm2 = sample(1:10000, Int(10000*0.2))
        p[rm1, 1] .= NaN; t[rm1, 1] .= NaN
        p[rm2, 2] .= NaN; t[rm2, 2] .= NaN
        # Compute the correlation
        MMasyn[i] = NUFFTcorrDKFGG(p, t)[1][1,2]
        HYasyn[i] = HYcorr(p[:,1], p[:,2], t[:,1], t[:,2])[1][1,2]
    end
    return MM, HY, MMasyn, HYasyn
end

# The function uses an OU process.
# Asynchrony is achieved by down-sampling by 20%.
function MeanRev(num)
    # Setup
    cor = range(-0.99, 0.99, length = num)
    MM = zeros(num, 1)
    HY = zeros(num, 1)
    MMasyn = zeros(num, 1)
    HYasyn = zeros(num, 1)
    # Parameters
    μ = [100, 100]
    σ1 = 0.1/86400
    σ2 = 0.2/86400
    θ = [0.035, 0.054]
    seed = 1:num
    # Run through each case
    @showprogress "Computing..." for i in 1:num
        # Simulate the GARCH
        Σ = [σ1 sqrt(σ1*σ2)*cor[i]; sqrt(σ1*σ2)*cor[i] σ2]
        p = OU(10000, μ, Σ, θ; seed = seed[i])
        t = [collect(1:1:10000.0) collect(1:1:10000.0)]
        # Compute the correlation
        MM[i] = NUFFTcorrDKFGG(p, t)[1][1,2]
        HY[i] = HYcorr(p[:,1], p[:,2], t[:,1], t[:,2])[1][1,2]
        # Downsample
        Random.seed!(seed[i])
        rm1 = sample(1:10000, Int(10000*0.2))
        Random.seed!(seed[i] + num)
        rm2 = sample(1:10000, Int(10000*0.2))
        p[rm1, 1] .= NaN; t[rm1, 1] .= NaN
        p[rm2, 2] .= NaN; t[rm2, 2] .= NaN
        # Compute the correlation
        MMasyn[i] = NUFFTcorrDKFGG(p, t)[1][1,2]
        HYasyn[i] = HYcorr(p[:,1], p[:,2], t[:,1], t[:,2])[1][1,2]
    end
    return MM, HY, MMasyn, HYasyn
end

## Compute the results
# Case with GARCH(1,1)
res1 = VolClust(21)

# Case with OU
res2 = MeanRev(21)

## Save the results
save("Computed Data/Reversion.jld", "res1", res1, "res2", res2)

## Load the results
res = load("Computed Data/Reversion.jld")
res1 = res["res1"]; res2 = res["res2"]

## Plot the results

p1 = plot(1:21, range(-0.99, 0.99, length = 21), legend = :topleft, label = "Induced", color = :black, linestyle = :dot, dpi = 300)
plot!(p1, 1:21, res1[1], label = "MM Syn", color = :blue, linestyle = :dash)
plot!(p1, 1:21, res1[2], label = "HY Syn", color = :red, linestyle = :dashdot)
plot!(p1, 1:21, res1[3], label = "MM Asyn", color = :purple, linestyle = :dash)
plot!(p1, 1:21, res1[4], label = "HY Asyn", color = :orange, linestyle = :dashdot)
ylabel!(p1, L"\textrm{Correlation }(\rho)")
xlabel!(p1, L"\textrm{Simulation}")
# savefig(p1, "Plots/RevGarch.svg")

p2 = plot(1:21, range(-0.99, 0.99, length = 21), legend = :topleft, label = "Induced", color = :black, linestyle = :dot, dpi = 300)
plot!(p2, 1:21, res2[1], label = "MM Syn", color = :blue, linestyle = :dash)
plot!(p2, 1:21, res2[2], label = "HY Syn", color = :red, linestyle = :dashdot)
plot!(p2, 1:21, res2[3], label = "MM Asyn", color = :purple, linestyle = :dash)
plot!(p2, 1:21, res2[4], label = "HY Asyn", color = :orange, linestyle = :dashdot)
ylabel!(p2, L"\textrm{Correlation }(\rho)")
xlabel!(p2, L"\textrm{Simulation}")
# savefig(p2, "Plots/RevOU.svg")
