## Author: Patrick Chang
# Script file to investigate asynchrony using the case
# of exponential inter-arrivals (Poisson sampling).
# Consider all the various models used before, and look
# at the impact of N in the Malliavin-Mancino estimator

## Preamble

using Plots, Random, ProgressMeter, StatsBase, JLD, LaTeXStrings, Distributions

cd("/Users/patrickchang1/PCRBTG-VT")

# Estimators
include("../Functions/Correlation Estimators/Dirichlet/NUFFTcorrDK-FGG.jl")
include("../Functions/Correlation Estimators/HY/HYcorr.jl")

# SDE
include("../Functions/SDEs/GBM.jl")
include("../Functions/SDEs/Merton Model.jl")
include("../Functions/SDEs/Variance Gamma.jl")
include("../Functions/SDEs/Garch.jl")
include("../Functions/SDEs/OU.jl")

## Supporting function to simulate random exponentials
function rexp(n, mean)
    t = -mean .* log.(rand(n))
end

## Function to streamline the process
# Here asynchrony is introduced using exponential inter-arrivals
# we sample the first price path with a mean of 30s and the second
# price path with a mean of 45s.
function asynchrony(process::Symbol)
    # Setup
    reps = 1000
    N = collect(20:10:200)
    MM = zeros(reps, length(N))
    HY = zeros(reps, 1)
    MMasyn = zeros(reps, length(N))
    HYasyn = zeros(reps, 1)
    # Initialise parameters depending on process
    if process == :GBM
        μ = [0.01/86400, 0.01/86400]
        σ1 = 0.1/86400; σ2 = 0.2/86400
    elseif process == :Merton
        μ = [0.01/86400, 0.01/86400]
        σ1 = 0.1/86400; σ2 = 0.2/86400
        a = [0,0]; b = [100/86400, 100/86400]
        λ = [0.2, 0.2]
    elseif process == :VG
        μ = [0.01/86400, 0.01/86400]
        σ1 = 0.1/86400; σ2 = 0.2/86400
        β = 1
    elseif process == :GARCH
        θ = [0.035, 0.054]
        λ = [0.296, 0.48]
        ω = [0.636, 0.476]
    elseif process == :OU
        μ = [100, 100]
        σ1 = 0.1/86400
        σ2 = 0.2/86400
        θ = [0.035, 0.054]
    else
        return println("Please input a valid process (check spelling).")
    end
    # Loop through the various replications
    @showprogress "Computing..." for i in 1:reps
        # Simulate the process
        if process == :GBM
            Σ = [σ1 sqrt(σ1*σ2)*0.35; sqrt(σ1*σ2)*0.35 σ2]
            p = GBM(10000, μ, Σ; seed = i)
        elseif process == :Merton
            Σ = [σ1 sqrt(σ1*σ2)*0.35; sqrt(σ1*σ2)*0.35 σ2]
            p = Merton(10000, μ, Σ, lambda, a, b; seed = i)
        elseif process == :VG
            # Simulate the VG
            Σ = [σ1 sqrt(σ1*σ2)*0.35; sqrt(σ1*σ2)*0.35 σ2]
            p = VG(10000, μ, Σ, β; seed = i)
        elseif process == :GARCH
            p = Garch(10000, θ, λ, ω, 0.35; seed = i)
        elseif process == :OU
            Σ = [σ1 sqrt(σ1*σ2)*0.35; sqrt(σ1*σ2)*0.35 σ2]
            p = OU(10000, μ, Σ, θ; seed = i)
        end
        # Synchronous times
        t = [collect(1:1:10000.0) collect(1:1:10000.0)]
        # Asynchronous samples
        # Sample times
        Random.seed!(i)
        t1 = [0; rexp(10000, 30)]
        t1 = cumsum(t1)
        t1 = filter((x) -> x < 10000, t1)
        Random.seed!(i+reps)
        t2 = [0; rexp(10000, 45)]
        t2 = cumsum(t2)
        t2 = filter((x) -> x < 10000, t2)
        # Sample the proces
        p1_asyn = p[Int.(floor.(t1)) .+ 1, 1]
        p2_asyn = p[Int.(floor.(t2)) .+ 1, 2]
        # Pad the vectors
        D = maximum([length(t1); length(t2)]) - minimum([length(t1); length(t2)])
        if length(t1) < length(t2)
            t1 = [t1; repeat([NaN], D)]
            p1_asyn = [p1_asyn; repeat([NaN], D)]
        else
            t2 = [t2; repeat([NaN], D)]
            p2_asyn = [p2_asyn; repeat([NaN], D)]
        end
        # Asynchronous times and prices
        t_asyn = [t1 t2]
        p_asyn = [p1_asyn p2_asyn]
        # Compute the HY results
        HY[i] = HYcorr(p[:,1], p[:,2], t[:,1], t[:,2])[1][1,2]
        HYasyn[i] = HYcorr(p_asyn[:,1], p_asyn[:,2], t_asyn[:,1], t_asyn[:,2])[1][1,2]
        # Compute the MM results
        for j in eachindex(N)
            MM[i,j] = NUFFTcorrDKFGG(p, t; N = N[j])[1][1,2]
            MMasyn[i,j] = NUFFTcorrDKFGG(p_asyn, t_asyn; N = N[j])[1][1,2]
        end
    end
    return MM, HY, MMasyn, HYasyn
end


## Compute the results
# Case with GBM
res1 = asynchrony(:GBM)

# Case with Merton
res2 = asynchrony(:Merton)

# Case with VG
res3 = asynchrony(:VG)

# Case with GARCH
res4 = asynchrony(:GARCH)

# Case with OU
res5 = asynchrony(:OU)

## Save the results
save("Computed Data/SamplingInterval.jld", "res1", res1, "res2", res2, "res3", res3, "res4", res4, "res5", res5)

## Load the results
res = load("Computed Data/SamplingInterval.jld")
res1 = res["res1"]; res2 = res["res2"]; res3 = res["res3"]; res4 = res["res4"]; res5 = res["res5"]

## Plot the results

reps = 1000
N = collect(20:10:200)
q = quantile.(TDist(reps-1), [0.975])
Nyq = 10000 / (2*30)

# GBM
err_MM = (q .* std(res1[1], dims = 1))
err_HY = (q .* std(res1[2]))
err_MMasyn = (q .* std(res1[3], dims = 1))
err_HYasyn = (q .* std(res1[4]))

p1 = plot(N, mean(res1[1], dims=1)', ribbon=err_MM', fillalpha=.3, legend = :topright, color = :blue, line=(1, [:dash]), label = L"\textrm{MM Syn}", marker=([:circle :d],4,2,stroke(4,:blue)), dpi = 300)
plot!(p1, N, mean(res1[3], dims=1)', ribbon=err_MMasyn', fillalpha=0.3, color = :red, line=(1, [:dash]), label = L"\textrm{MM Asyn}", marker=([:x :d],4,2,stroke(4,:red)))
hline!(p1, N, [mean(res1[2])], ribbon=err_HY, fillalpha=.3, color = :green, line=(1, [:dashdot]), label = L"\textrm{HY Syn}", marker=([:circle :d],4,2,stroke(2,:green)))
hline!(p1, N, [mean(res1[4])], ribbon=err_HYasyn, fillalpha=0.5, color = :orange, line=(1, [:dashdot]), label = L"\textrm{HY Asyn}", marker=([:x :d],4,2,stroke(2,:orange)))
hline!(p1, [.35], color = :black, line=(2, [:dot]), label = L"\textrm{Induced } \rho")
vline!(p1, [Nyq], line=(2, [:dot]), color = :black)
annotate!(p1, 160, 0.03, text(L"\textrm{Average Nyquist cutoff}", :black, :right, 8))
xlabel!(p1, L"\textrm{# of Fourier modes } (N)")
ylabel!(p1, L"\textrm{Correlation } (\rho)")
# savefig(p1, "Plots/SIGBM.svg")

# Merton
err_MM = (q .* std(res2[1], dims = 1))
err_HY = (q .* std(res2[2]))
err_MMasyn = (q .* std(res2[3], dims = 1))
err_HYasyn = (q .* std(res2[4]))

p2 = plot(N, mean(res2[1], dims=1)', ribbon=err_MM', fillalpha=.3, legend = :topright, color = :blue, line=(1, [:dash]), label = L"\textrm{MM Syn}", marker=([:circle :d],4,2,stroke(4,:blue)), dpi = 300)
plot!(p2, N, mean(res2[3], dims=1)', ribbon=err_MMasyn', fillalpha=0.3, color = :red, line=(1, [:dash]), label = L"\textrm{MM Asyn}", marker=([:x :d],4,2,stroke(4,:red)))
hline!(p2, N, [mean(res2[2])], ribbon=err_HY, fillalpha=.3, color = :green, line=(1, [:dashdot]), label = L"\textrm{HY Syn}", marker=([:circle :d],4,2,stroke(2,:green)))
hline!(p2, N, [mean(res2[4])], ribbon=err_HYasyn, fillalpha=0.5, color = :orange, line=(1, [:dashdot]), label = L"\textrm{HY Asyn}", marker=([:x :d],4,2,stroke(2,:orange)))
hline!(p2, [.35], color = :black, line=(2, [:dot]), label = L"\textrm{Induced } \rho")
vline!(p2, [Nyq], line=(2, [:dot]), color = :black)
annotate!(p2, 160, 0.02, text(L"\textrm{Average Nyquist cutoff}", :black, :right, 8))
xlabel!(p2, L"\textrm{# of Fourier modes } (N)")
ylabel!(p2, L"\textrm{Correlation } (\rho)")
# savefig(p2, "Plots/SIMerton.svg")

# VG
err_MM = (q .* std(res3[1], dims = 1))
err_HY = (q .* std(res3[2]))
err_MMasyn = (q .* std(res3[3], dims = 1))
err_HYasyn = (q .* std(res3[4]))

p3 = plot(N, mean(res3[1], dims=1)', ribbon=err_MM', fillalpha=.3, legend = :topright, color = :blue, line=(1, [:dash]), label = L"\textrm{MM Syn}", marker=([:circle :d],4,2,stroke(4,:blue)), dpi = 300)
plot!(p3, N, mean(res3[3], dims=1)', ribbon=err_MMasyn', fillalpha=0.3, color = :red, line=(1, [:dash]), label = L"\textrm{MM Asyn}", marker=([:x :d],4,2,stroke(4,:red)))
hline!(p3, N, [mean(res3[2])], ribbon=err_HY, fillalpha=.3, color = :green, line=(1, [:dashdot]), label = L"\textrm{HY Syn}", marker=([:circle :d],4,2,stroke(2,:green)))
hline!(p3, N, [mean(res3[4])], ribbon=err_HYasyn, fillalpha=0.5, color = :orange, line=(1, [:dashdot]), label = L"\textrm{HY Asyn}", marker=([:x :d],4,2,stroke(2,:orange)))
hline!(p3, [.35], color = :black, line=(2, [:dot]), label = L"\textrm{Induced } \rho")
vline!(p3, [Nyq], line=(2, [:dot]), color = :black)
annotate!(p3, 160, 0.032, text(L"\textrm{Average Nyquist cutoff}", :black, :right, 8))
xlabel!(p3, L"\textrm{# of Fourier modes } (N)")
ylabel!(p3, L"\textrm{Correlation } (\rho)")
# savefig(p3, "Plots/SIVG.svg")

# GARCH
err_MM = (q .* std(res4[1], dims = 1))
err_HY = (q .* std(res4[2]))
err_MMasyn = (q .* std(res4[3], dims = 1))
err_HYasyn = (q .* std(res4[4]))

p4 = plot(N, mean(res4[1], dims=1)', ribbon=err_MM', fillalpha=.3, legend = :topright, color = :blue, line=(1, [:dash]), label = L"\textrm{MM Syn}", marker=([:circle :d],4,2,stroke(4,:blue)), dpi = 300)
plot!(p4, N, mean(res4[3], dims=1)', ribbon=err_MMasyn', fillalpha=0.3, color = :red, line=(1, [:dash]), label = L"\textrm{MM Asyn}", marker=([:x :d],4,2,stroke(4,:red)))
hline!(p4, N, [mean(res4[2])], ribbon=err_HY, fillalpha=.3, color = :green, line=(1, [:dashdot]), label = L"\textrm{HY Syn}", marker=([:circle :d],4,2,stroke(2,:green)))
hline!(p4, N, [mean(res4[4])], ribbon=err_HYasyn, fillalpha=0.5, color = :orange, line=(1, [:dashdot]), label = L"\textrm{HY Asyn}", marker=([:x :d],4,2,stroke(2,:orange)))
hline!(p4, [.35], color = :black, line=(2, [:dot]), label = L"\textrm{Induced } \rho")
vline!(p4, [Nyq], line=(2, [:dot]), color = :black)
annotate!(p4, 160, 0.04, text(L"\textrm{Average Nyquist cutoff}", :black, :right, 8))
xlabel!(p4, L"\textrm{# of Fourier modes } (N)")
ylabel!(p4, L"\textrm{Correlation } (\rho)")
# savefig(p4, "Plots/SIGARCH.svg")

# OU
err_MM = (q .* std(res5[1], dims = 1))
err_HY = (q .* std(res5[2]))
err_MMasyn = (q .* std(res5[3], dims = 1))
err_HYasyn = (q .* std(res5[4]))

p5 = plot(N, mean(res5[1], dims=1)', ribbon=err_MM', fillalpha=.3, legend = :topright, color = :blue, line=(1, [:dash]), label = L"\textrm{MM Syn}", marker=([:circle :d],4,2,stroke(4,:blue)), dpi = 300)
plot!(p5, N, mean(res5[3], dims=1)', ribbon=err_MMasyn', fillalpha=0.3, color = :red, line=(1, [:dash]), label = L"\textrm{MM Asyn}", marker=([:x :d],4,2,stroke(4,:red)))
hline!(p5, N, [mean(res5[2])], ribbon=err_HY, fillalpha=.3, color = :green, line=(1, [:dashdot]), label = L"\textrm{HY Syn}", marker=([:circle :d],4,2,stroke(2,:green)))
hline!(p5, N, [mean(res5[4])], ribbon=err_HYasyn, fillalpha=0.5, color = :orange, line=(1, [:dashdot]), label = L"\textrm{HY Asyn}", marker=([:x :d],4,2,stroke(2,:orange)))
hline!(p5, [.35], color = :black, line=(2, [:dot]), label = L"\textrm{Induced } \rho")
vline!(p5, [Nyq], line=(2, [:dot]), color = :black)
annotate!(p5, 160, 0, text(L"\textrm{Average Nyquist cutoff}", :black, :right, 8))
xlabel!(p5, L"\textrm{# of Fourier modes } (N)")
ylabel!(p5, L"\textrm{Correlation } (\rho)")
# savefig(p5, "Plots/SIOU.svg")
