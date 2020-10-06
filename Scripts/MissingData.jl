## Author: Patrick Chang
# Script file to investigate the effect of asynchrony
# introduced using the missing data approach

## Preamble

using Plots, Random, ProgressMeter, StatsBase, JLD, LaTeXStrings

cd("/Users/patrickchang1/PCRBTG-VT")

# Estimators
include("../Functions/Correlation Estimators/Dirichlet/NUFFTcorrDK-FGG.jl")
include("../Functions/Correlation Estimators/HY/HYcorr.jl")

# SDE
include("../Functions/SDEs/GBM.jl")

## Function to streamline the process
function MissingData(num, percentage)
    # Setup
    cor = range(-0.99, 0.99, length = num)
    MM = zeros(num, 1)
    HY = zeros(num, 1)
    # Parameters
    μ = [0.01/86400, 0.01/86400]
    σ1 = 0.1/86400; σ2 = 0.2/86400
    seed = 1:num
    # Run through each case
    @showprogress "Computing..." for i in 1:num
        # Simulate the GBM
        Σ = [σ1 sqrt(σ1*σ2)*cor[i]; sqrt(σ1*σ2)*cor[i] σ2]
        p = GBM(10000, μ, Σ; seed = seed[i])
        t = [collect(1:1:10000.0) collect(1:1:10000.0)]
        # Downsample
        Random.seed!(seed[i])
        rm1 = sample(1:10000, Int(10000*percentage))
        Random.seed!(seed[i] + num)
        rm2 = sample(1:10000, Int(10000*percentage))
        p[rm1, 1] .= NaN; t[rm1, 1] .= NaN
        p[rm2, 2] .= NaN; t[rm2, 2] .= NaN
        # Compute the correlation
        MM[i] = NUFFTcorrDKFGG(p, t)[1][1,2]
        HY[i] = HYcorr(p[:,1], p[:,2], t[:,1], t[:,2])[1][1,2]
    end
    return MM, HY
end


## Compute the results
# Case with no missing data
res1 = MissingData(21, 0)

# Case with 10% missing data
res2 = MissingData(21, 0.1)

# Case with 20% missing data
res3 = MissingData(21, 0.2)

# Case with 40% missing data
res4 = MissingData(21, 0.4)

## Save the results
save("Computed Data/MissingData.jld", "res1", res1, "res2", res2, "res3", res3, "res4", res4)

## Load the results
res = load("Computed Data/MissingData.jld")
res1 = res["res1"]; res2 = res["res2"]; res3 = res["res3"]; res4 = res["res4"]

## Plot the results

p1 = plot(1:21, range(-0.99, 0.99, length = 21), legend = :topleft, label = "Induced", color = :black, linestyle = :dot, dpi = 300)
plot!(p1, 1:21, res1[1], label = "MM", color = :blue, linestyle = :dash)
plot!(p1, 1:21, res1[2], label = "HY", color = :red, linestyle = :dashdot)
ylabel!(p1, L"\textrm{Correlation }(\rho)")
xlabel!(p1, L"\textrm{Simulation}")
# savefig(p1, "Plots/MissingData0.svg")

p2 = plot(1:21, range(-0.99, 0.99, length = 21), legend = :topleft, label = "Induced", color = :black, linestyle = :dot, dpi = 300)
plot!(p2, 1:21, res2[1], label = "MM", color = :blue, linestyle = :dash)
plot!(p2, 1:21, res2[2], label = "HY", color = :red, linestyle = :dashdot)
ylabel!(p2, L"\textrm{Correlation }(\rho)")
xlabel!(p2, L"\textrm{Simulation}")
# savefig(p2, "Plots/MissingData1.svg")

p3 = plot(1:21, range(-0.99, 0.99, length = 21), legend = :topleft, label = "Induced", color = :black, linestyle = :dot, dpi = 300)
plot!(p3, 1:21, res3[1], label = "MM", color = :blue, linestyle = :dash)
plot!(p3, 1:21, res3[2], label = "HY", color = :red, linestyle = :dashdot)
ylabel!(p3, L"\textrm{Correlation }(\rho)")
xlabel!(p3, L"\textrm{Simulation}")
# savefig(p3, "Plots/MissingData2.svg")

p4 = plot(1:21, range(-0.99, 0.99, length = 21), legend = :topleft, label = "Induced", color = :black, linestyle = :dot, dpi = 300)
plot!(p4, 1:21, res4[1], label = "MM", color = :blue, linestyle = :dash)
plot!(p4, 1:21, res4[2], label = "HY", color = :red, linestyle = :dashdot)
ylabel!(p4, L"\textrm{Correlation }(\rho)")
xlabel!(p4, L"\textrm{Simulation}")
# savefig(p4, "Plots/MissingData4.svg")
