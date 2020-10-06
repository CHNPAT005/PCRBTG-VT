## Author: Patrick Chang
# Function to simulate a 2 dimensional Bates model
# using the first order Euler discretization

using Distributions, Random, LinearAlgebra, Plots, LaTeXStrings

#---------------------------------------------------------------------------
# ρ_σ = the correlation between the volatility process
# ρ = the correlation between the volatility process and the price process
# λ = the vector of rate of jumps on the volatility and price process
# ν_p = the vector of variance of jump size of price process
# ν_σ = the vector of variance of jump size of volatility process

# starting price is at 100 and volatility is at 0.04

# parameters are preset using the choices from Cuchiero and Teichmann 2015
function Bates_CT(n, ρ = [-0.3; -0.5], M = [-1.6 -0.2; -0.4 -1], α = [0.0725 0.06; 0.06 0.1325], λ_Y = 100, λ_X = 10, μ = [-0.005; -0.003], ν = [0.015; 0.02], θ = 0.05; kwargs...)
    kwargs = Dict(kwargs)

    if haskey(kwargs, :SP)
        SP = kwargs[:SP]
    else
        SP = 100
    end

    if haskey(kwargs, :seed)
        seed = kwargs[:seed]
    else
        seed = 1
    end

    if haskey(kwargs, :SV)
        SV = kwargs[:SV]
    else
        SV = sqrt(0.04)
    end

    local dt::Float64
    if haskey(kwargs, :dt)
        dt = kwargs[:dt]
    else
        dt = 21600
    end

    X_minus = zeros(ComplexF64, 2, 2, n)
    X_minus[:,:,1] = [0.09 -0.036; -0.036 0.09]
    X_plus = zeros(ComplexF64, 2, 2, n)
    X_plus[:,:,1] = [0.09 -0.036; -0.036 0.09]
    Y = zeros(ComplexF64, n, 2)
    Y[1,:] = log.([100 100])

    Σ = sqrt(α)
    b = 3.5*α

    Random.seed!(seed)
    W = randn(2, n)
    Random.seed!(2*seed)
    B = randn(2, 2, n)
    Random.seed!(3*seed)
    Z_jump = randn(2, n)

    Random.seed!(seed)
    N_X = rand(Poisson(λ_X/dt), n)
    Random.seed!(2*seed)
    N_Y1 = rand(Poisson(λ_Y/dt), n)
    Random.seed!(3*seed)
    N_Y2 = rand(Poisson(λ_Y/dt), n)
    N_Y = [N_Y1 N_Y2]

    for i in 2:n
        X_minus[:,:,i] = X_plus[:,:,i-1] + (b + M*X_plus[:,:,i-1] + X_plus[:,:,i-1]*M')/dt + sqrt(X_plus[:,:,i-1])*B[:,:,i]*Σ/sqrt(dt) + Σ*B[:,:,i]'*sqrt(X_plus[:,:,i-1])/sqrt(dt)
        if N_X[i] == 0
            X_plus[:,:,i] = X_minus[:,:,i]
        else
            X_plus[:,:,i] = [X_minus[:,:,i][1,1]+sum(rand(Exponential(θ), N_X[i])) X_minus[:,:,i][1,2]; X_minus[:,:,i][2,1] X_minus[:,:,i][2,2]]
        end
    end

    for i in 2:n
        Z = sqrt(1-ρ'ρ)*W[:,i] + B[:,:,i]*ρ
        bs = -0.5*diag(X_plus[:,:,i]) - λ_Y.*(exp.(μ - 0.5*ν.^2) .- 1)
        M = μ.*N_Y[i,:] + ν.*sqrt.(N_Y[i,:]).*Z_jump[:,i]
        Y[i,:] = Y[i-1,:] + bs/dt + sqrt(X_minus[:,:,i])*Z/sqrt(dt) + M
    end

    σ_1 = zeros(ComplexF64, n, 1)
    σ_2 = zeros(ComplexF64, n, 1)
    σ_12 = zeros(ComplexF64, n, 1)

    for i in 1:n
        σ_1[i] = X_plus[:,:,i][1,1]
        σ_2[i] = X_plus[:,:,i][2,2]
        σ_12[i] = X_plus[:,:,i][1,2]
    end

    return real(exp.(Y)), real(σ_1), real(σ_2), real(σ_12)
end

# theta = [0.035 0.054]
# lambda  = [0.296 0.48]
# w = [0.636 0.476]
# aa = GARCH_Reno(21600, theta, lambda, w, 0.35, dt = 21600)
#
# aa = Heston_CT(21600, seed = 2)#, dt = 86400)
# # aa = Heston2D(511*250, seed = 10, dt = 511)
# t = collect(1:1:21600)
# a = NUFFTcorrDKFGG(aa[1], [t t])
# a[1]
#
# p2 = plot(t, real(aa[2]), label = L"\sigma_1")
# plot!(p2, t, real(aa[3]), label = L"\sigma_2")
# plot!(p2, t, real(aa[4]), label = L"\sigma_{12}")
#
#
#
# #---------------------------------------------------------------------------
#
# aa = Heston2D(21600, seed = 1)
# aa = Bates2D(21600, seed = 100)
#
# aa = Heston_CT(21600, seed = 2)
#
# y1 = log.(real(aa[1][:,1]))
# y2 = log.(real(aa[1][:,2]))
#
# dj_bates1 = diff(y1, dims = 1)
# dj_bates2 = diff(y2, dims = 1)
# dj_bates12 = diff(y1+y2, dims = 1)
#
# N = 128
#
# M = N*2+1
#
# # T = collect(1:1:M)
# # T = (T .- minimum(T)) .* (2*pi / (maximum(T) - minimum(T)))
# t = collect(1:1:21600)
# t = (t .- minimum(t)) .* (2*pi / (maximum(t) - minimum(t)))
# t = t[1:end-1]
# T = t[Int.(floor.(collect(1:21600/M:21600)))]
#
#
# rho_bates1 = test(T, dj_bates1, t, N)
# rho_bates2 = test(T, dj_bates2, t, N)
# rho_bates12 = test(T, dj_bates12, t, N)
#
# t = collect(1:1:21600)
# t = (t .- minimum(t)) .* (1 / (maximum(t) - minimum(t)))
# t = t[1:end-1]
#
# rho_bates1 = hat_ρ(dj_bates1, N, 21600-1, t)
# rho_bates2 = hat_ρ(dj_bates2, N, 21600-1, t)
# rho_bates12 = hat_ρ(dj_bates12, N, 21600-1, t)
#
#
# sighat_bates1 = -2 .* log.(real(rho_bates1))
# sighat_bates2 = -2 .* log.(real(rho_bates2))
# sighat_bates12 = 0.5 .* ( (-2 .* log.(real(rho_bates12))) - sighat_bates1 - sighat_bates2 )
#
#
# t = collect(1:1:21600)
# t = (t .- minimum(t)) .* (1 / (maximum(t) - minimum(t)))
#
# b1 = MM_JR(aa[1], N)[1]
# b1 = LRV(aa[1], M)[1]
# # σ_1
# a = T./(2*pi)
# b1 = sighat_bates1
# p1 = plot(t, real(aa[2])[1:end-1], label = L"\textrm{Simulated } \sigma_{11}^2", line=(0.5, [:solid]))
# plot!(p1, t, b1, label = L"\textrm{Estimated } \sigma_{11}^2")
# # p1 = plot(a, real(aa[2][Int.(floor.(collect(1:21600/M:21600)))]), label = L"\textrm{Simulated } \sigma_{11}^2")
# # plot!(p1, a, b1, label = L"\textrm{Estimated } \sigma_{11}^2")
# xlabel!(p1, "Time")
# ylabel!(p1, L"\sigma_{11}^2(t)")
#
# # savefig(p1, "Plots/HestonVol11.pdf")
#
# # σ_2
# a = T./(2*pi)
# b = sighat_bates2
# p2 = plot(t, real(aa[3]), label = L"\textrm{Simulated } \sigma_{22}^2", line=(0.5, [:solid]))
# plot!(p2, a[1:end-1], b[1:end-1], label = L"\textrm{Estimated } \sigma_{22}^2")
# # p2 = plot(a, real(aa[3][Int.(floor.(collect(1:21600/M:21600)))]), label = L"\textrm{Simulated } \sigma_{11}^2")
# # plot!(p2, a, b, label = L"\textrm{Estimated } \sigma_{11}^2")
# xlabel!(p2, "Time")
# ylabel!(p2, L"\sigma_2^{22}(t)")
#
# # savefig(p2, "Plots/HestonVol22.pdf")
#
# # σ_12
# a = T./(2*pi)
# b2 = sighat_bates12
# p3 = plot(t, aa[4], label = L"\textrm{Simulated } \sigma_{12}^2", line=(0.5, [:solid]))
# plot!(p3, a[1:end-1], b2[1:end-1], label = L"\textrm{Estimated } \sigma_{12}^2")
# # p3 = plot(a, real(aa[4][Int.(floor.(collect(1:21600/M:21600)))]), label = L"\textrm{Simulated } \sigma_{11}^2")
# # plot!(p3, a, b2, label = L"\textrm{Estimated } \sigma_{11}^2")
# xlabel!(p3, "Time")
# ylabel!(p3, L"\sigma_{12}^2(t)")
#
# # savefig(p3, "Plots/HestonVol12.pdf")
#
# p4 = plot(t, aa[4] ./ (sqrt.(aa[2]) .* sqrt.(aa[3])), label = L"\textrm{Simulated } \rho")
# plot!(p4, a[1:end-1], b2[1:end-1] ./ (sqrt.(b1[1:end-1]) .* sqrt.(b[1:end-1])), label = L"\textrm{Estimated } \rho")
#
#
#
# p1 = plot(t, real(aa[1][:,1]), label = L"P_1")
# plot!(p1, t, real(aa[1][:,2]), label = L"P_2")
#
# p2 = plot(t, aa[2], label = L"\sigma_{11}")
# plot!(p2, t, aa[3], label = L"\sigma_{22}")
# plot!(p2, t, aa[4], label = L"\sigma_{12}")
# xlabel!(p2, "Time")
# ylabel!(p2, L"\sigma_{ij}^2(t)")
#
# # savefig(p2, "Plots/HestonVol.pdf")
#
# # p3 = plot(t, aa[4], label = L"\sigma_{12}")
# #---------------------------------------------------------------------------
#
# p1 = plot(t, aa[2], label = "Simulated", color = :blue, line=(0.5, [:solid]))
# plot!(p1, t, aa[3], label = "", color = :blue, line=(0.5, [:solid]))
# plot!(p1, t, aa[4], label = "", color = :blue, line=(0.5, [:solid]))
# plot!(p1, a[1:end-1], b1[1:end-1], label = "Estimated", color = :orange)
# plot!(p1, a[1:end-1], b[1:end-1], label = "", color = :orange)
# plot!(p1, a[1:end-1], b2[1:end-1], label = "", color = :orange)
#
#
# #---------------------------------------------------------------------------
# mu = [0.01/21600, 0.01/21600]
# sigma = [0.1/21600 sqrt(0.1/21600)*0.35*sqrt(0.2/21600);
#         sqrt(0.1/21600)*0.35*sqrt(0.2/21600) 0.2/21600]
#
#
# t = collect(1:1:21600)
# aa = Heston_CT(21600, seed = 2)
# aa = GBM(21600, mu, sigma)
# # b = MM_JR(aa[1], 128, 1000)
# b = MM_inst(aa, [t t], 100)
#
#
# t = collect(1:1:21600)
# t = (t .- minimum(t)) .* (1 / (maximum(t) - minimum(t)))
#
# tt = collect(0:1/100:1-(1/100))
# p1 = plot(t, aa[2], label = L"\textrm{Simulated } \sigma_{11}^2", color = :blue, line=(0.5, [:solid]))
# plot!(p1, tt, b[1], label = L"\textrm{Estimated } \sigma_{11}^2", color = :orange, line=(1, [:solid]))
#
# p1 = plot(t, aa[3], label = L"\textrm{Simulated } \sigma_{22}^2", color = :blue, line=(0.5, [:solid]))
# plot!(p1, tt, b[2], label = L"\textrm{Estimated } \sigma_{22}^2", color = :orange, line=(1, [:solid]))
#
# p1 = plot(t, aa[4], label = L"\textrm{Simulated } \sigma_{12}^2", color = :blue, line=(0.5, [:solid]))
# plot!(p1, tt, b[3], label = L"\textrm{Estimated } \sigma_{12}^2", color = :orange, line=(1, [:solid]))
#
#
# p1 = plot(tt, b[1], label = L"\textrm{Estimated } \sigma_{11}^2", color = :orange, line=(1, [:solid]))
# hline!(p1, [0.1])
# p1 = plot(tt, b[2], label = L"\textrm{Estimated } \sigma_{11}^2", color = :orange, line=(1, [:solid]))
# hline!(p1, [0.2])
# p1 = plot(tt, b[3], label = L"\textrm{Estimated } \sigma_{11}^2", color = :orange, line=(1, [:solid]))
# hline!(p1, [sqrt(0.1/1)*0.35*sqrt(0.2/1)])
