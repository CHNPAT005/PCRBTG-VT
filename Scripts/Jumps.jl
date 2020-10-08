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
include("../Functions/SDEs/Merton Model.jl")
include("../Functions/SDEs/Variance Gamma.jl")

## Function to streamline the process
# The function uses a Merten model with λ
# as a varying input. Here for the asynchronous Case
# we down-sample using 20%.
function Jumps(num, λ)
    # Setup
    cor = range(-0.99, 0.99, length = num)
    MM = zeros(num, 1)
    HY = zeros(num, 1)
    MMasyn = zeros(num, 1)
    HYasyn = zeros(num, 1)
    # Parameters
    μ = [0.01/86400, 0.01/86400]
    σ1 = 0.1/86400; σ2 = 0.2/86400
    a = [0,0]
    b = [100/86400, 100/86400]
    lambda = [λ, λ]
    seed = 1:num
    # Run through each case
    @showprogress "Computing..." for i in 1:num
        # Simulate the GBM
        Σ = [σ1 sqrt(σ1*σ2)*cor[i]; sqrt(σ1*σ2)*cor[i] σ2]
        p = Merton(10000, μ, Σ, lambda, a, b; seed = seed[i])
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

# The function uses a Variance Gamma model.
# Again the asynchronous cases we down-sample using 20%.
function PureJumps(num)
    # Setup
    cor = range(-0.99, 0.99, length = num)
    MM = zeros(num, 1)
    HY = zeros(num, 1)
    MMasyn = zeros(num, 1)
    HYasyn = zeros(num, 1)
    # Parameters
    μ = [0.01/86400, 0.01/86400]
    σ1 = 0.1/86400; σ2 = 0.2/86400
    β = 1
    seed = 1:num
    # Run through each case
    @showprogress "Computing..." for i in 1:num
        # Simulate the GBM
        Σ = [σ1 sqrt(σ1*σ2)*cor[i]; sqrt(σ1*σ2)*cor[i] σ2]
        p = VG(10000, μ, Σ, β; seed = seed[i])
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
# Case with λ = 0
res1 = Jumps(21, 0)

# Case with λ = 0.2
res2 = Jumps(21, 0.2)

# Case with λ = 0.5
res3 = Jumps(21, 0.5)

# Case with pure jumps (Vairance Gamma)
res4 = PureJumps(21)

## Save the results
save("Computed Data/Jumps.jld", "res1", res1, "res2", res2, "res3", res3, "res4", res4)

## Load the results
res = load("Computed Data/Jumps.jld")
res1 = res["res1"]; res2 = res["res2"]; res3 = res["res3"]; res4 = res["res4"]

## Plot the results

p1 = plot(1:21, range(-0.99, 0.99, length = 21), legend = :topleft, label = "Induced", color = :black, linestyle = :dot, dpi = 300)
plot!(p1, 1:21, res1[1], label = "MM Syn", color = :blue, linestyle = :dash)
plot!(p1, 1:21, res1[2], label = "HY Syn", color = :red, linestyle = :dashdot)
plot!(p1, 1:21, res1[3], label = "MM Asyn", color = :purple, linestyle = :dash)
plot!(p1, 1:21, res1[4], label = "HY Asyn", color = :orange, linestyle = :dashdot)
ylabel!(p1, L"\textrm{Correlation }(\rho)")
xlabel!(p1, L"\textrm{Simulation}")
# savefig(p1, "Plots/Jumps0.svg")

p2 = plot(1:21, range(-0.99, 0.99, length = 21), legend = :topleft, label = "Induced", color = :black, linestyle = :dot, dpi = 300)
plot!(p2, 1:21, res2[1], label = "MM Syn", color = :blue, linestyle = :dash)
plot!(p2, 1:21, res2[2], label = "HY Syn", color = :red, linestyle = :dashdot)
plot!(p2, 1:21, res2[3], label = "MM Asyn", color = :purple, linestyle = :dash)
plot!(p2, 1:21, res2[4], label = "HY Asyn", color = :orange, linestyle = :dashdot)
ylabel!(p2, L"\textrm{Correlation }(\rho)")
xlabel!(p2, L"\textrm{Simulation}")
# savefig(p2, "Plots/Jumps2.svg")

p3 = plot(1:21, range(-0.99, 0.99, length = 21), legend = :topleft, label = "Induced", color = :black, linestyle = :dot, dpi = 300)
plot!(p3, 1:21, res3[1], label = "MM Syn", color = :blue, linestyle = :dash)
plot!(p3, 1:21, res3[2], label = "HY Syn", color = :red, linestyle = :dashdot)
plot!(p3, 1:21, res3[3], label = "MM Asyn", color = :purple, linestyle = :dash)
plot!(p3, 1:21, res3[4], label = "HY Asyn", color = :orange, linestyle = :dashdot)
ylabel!(p3, L"\textrm{Correlation }(\rho)")
xlabel!(p3, L"\textrm{Simulation}")
# savefig(p3, "Plots/Jumps5.svg")

p4 = plot(1:21, range(-0.99, 0.99, length = 21), legend = :topleft, label = "Induced", color = :black, linestyle = :dot, dpi = 300)
plot!(p4, 1:21, res4[1], label = "MM Syn", color = :blue, linestyle = :dash)
plot!(p4, 1:21, res4[2], label = "HY Syn", color = :red, linestyle = :dashdot)
plot!(p4, 1:21, res4[3], label = "MM Asyn", color = :purple, linestyle = :dash)
plot!(p4, 1:21, res4[4], label = "HY Asyn", color = :orange, linestyle = :dashdot)
ylabel!(p4, L"\textrm{Correlation }(\rho)")
xlabel!(p4, L"\textrm{Simulation}")
# savefig(p4, "Plots/JumpsPure.svg")
