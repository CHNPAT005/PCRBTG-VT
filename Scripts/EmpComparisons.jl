## Author: Patrick Chang
# Script file to investigate the Epps effect in calendar time
# for 4 banking stocks on the JSE.

## Preamble

using CSV, Plots, Random, ProgressMeter, StatsBase, JLD, LaTeXStrings, DataFrames, Dates, Distributions, Pipe, ColorSchemes

cd("/Users/patrickchang1/PCRBTG-VT")

## Body of paper
#--------------------------------------------------------------------------------
## Read in the data to get the ticker names
prices = CSV.read("Real Data/JSE_prices_2019-06-24_2019-06-28.csv")
# Remove extra column
prices = prices[:,2:end]
# Pull out banking stocks
prices = prices[:,[6:8; 10]]
# Get the names of tickers
tickers = names(prices)

## Load results
CT = load("Computed Data/EmpCT.jld")
CT = CT["res"]

ET = load("Computed Data/EmpET.jld")
ET = ET["res"]

VT = load("Computed Data/EmpVT.jld")
VT = VT["res"]

## Extract label pairs
pairnames = Matrix{Union{Nothing, String}}(nothing, 1, 6)
inds = 1
for i in 1:(length(tickers)-1)
    for j in (i+1):length(tickers)
        # push!(pairnames, "$(tickers[i])"*"/"*"$(tickers[j])")
        pairnames[inds] = "$(tickers[i])"*"/"*"$(tickers[j])"
        inds += 1
    end
end

## Plot the results
dt = collect(1:1:300)

# SBK/FSR pair (at index number 3)
plot(dt, CT[1][:,3], legend = :right, color = :blue, line=(1, [:dash]), label = L"\textrm{CT MM}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)))
plot!(dt, CT[2][:,3], color = :red, line=(1, [:dash]), label = L"\textrm{CT RV}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)))
plot!(dt, ET[1][:,3], color = :green, line=(1, [:dash]), label = L"\textrm{ET MM}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)))
plot!(dt, ET[2][:,3], color = :orange, line=(1, [:dash]), label = L"\textrm{ET RV}", marker = ([:x :d], 3, 0.8, stroke(3, :red)))
plot!(dt, VT[:,3], color = :pink, line=(1, [:dash]), label = L"\textrm{VT RV}", marker = ([:+ :d], 3, 0.8, stroke(3, :red)))
hline!([CT[3][3]], color = :brown, line=(1, [:dash]), label = L"\textrm{CT HY}")
hline!([ET[3][3]], color = :lightblue, line=(2, [:dashdot]), label = L"\textrm{ET HY}")
xlabel!(L"\textrm{Sampling interval}")
ylabel!(L"\tilde{\rho}_{\Delta t}^{ij}")
# savefig("Plots/EmpCompSBKFSR.svg")


# NED/ABG pair (at index number 4)
plot(dt, CT[1][:,4], legend = :right, color = :blue, line=(1, [:dash]), label = L"\textrm{CT MM}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)))
plot!(dt, CT[2][:,4], color = :red, line=(1, [:dash]), label = L"\textrm{CT RV}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)))
plot!(dt, ET[1][:,4], color = :green, line=(1, [:dash]), label = L"\textrm{ET MM}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)))
plot!(dt, ET[2][:,4], color = :orange, line=(1, [:dash]), label = L"\textrm{ET RV}", marker = ([:x :d], 3, 0.8, stroke(3, :red)))
plot!(dt, VT[:,4], color = :pink, line=(1, [:dash]), label = L"\textrm{VT RV}", marker = ([:+ :d], 3, 0.8, stroke(3, :red)))
hline!([CT[3][4]], color = :brown, line=(1, [:dash]), label = L"\textrm{CT HY}")
hline!([ET[3][4]], color = :lightblue, line=(2, [:dashdot]), label = L"\textrm{ET HY}")
xlabel!(L"\textrm{Sampling interval}")
ylabel!(L"\tilde{\rho}_{\Delta t}^{ij}")
# savefig("Plots/EmpCompNEDABG.svg")
#--------------------------------------------------------------------------------


## Appendix of paper
#--------------------------------------------------------------------------------
## Read in the data to get the ticker names
prices = CSV.read("Real Data/JSE_prices_2019-06-24_2019-06-28.csv", DataFrame)
# Remove extra column
prices = prices[:,2:end]
# Pull out banking stocks
prices = prices[:,[1;3;6;7]]
# Get the names of tickers
tickers = names(prices)

## Load results
CTMM = load("Computed Data/EmpAppCTMM.jld")
CTMM = CTMM["res"]
CTRV = load("Computed Data/EmpAppCTRV.jld")
CTRV = CTRV["res"]
CTHY = load("Computed Data/EmpAppCTHY.jld")
CTHY = CTHY["res"]

ETMM = load("Computed Data/EmpAppETMM.jld")
ETMM = ETMM["res"]
ETRV = load("Computed Data/EmpAppETRV.jld")
ETRV = ETRV["res"]
ETHY = load("Computed Data/EmpAppETHY.jld")
ETHY = ETHY["res"]

VT = load("Computed Data/EmpAppVT.jld")
VT = VT["res"]

## Extract label pairs
pairnames = Matrix{Union{Nothing, String}}(nothing, 1, 6)
inds = 1
for i in 1:(length(tickers)-1)
    for j in (i+1):length(tickers)
        # push!(pairnames, "$(tickers[i])"*"/"*"$(tickers[j])")
        pairnames[inds] = "$(tickers[i])"*"/"*"$(tickers[j])"
        inds += 1
    end
end

## Plot the results
dt = collect(1:1:300)

# BTI/NED pair (at index number 3)
plot(dt, CTMM[:,3], legend = :bottomleft, color = :blue, line=(1, [:dash]), label = L"\textrm{CT MM}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)), ylims = (-0.3, 0.3))
plot!(dt, CTRV[:,3], color = :red, line=(1, [:dash]), label = L"\textrm{CT RV}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)))
plot!(dt, ETMM[:,3], color = :green, line=(1, [:dash]), label = L"\textrm{ET MM}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)))
plot!(dt, ETRV[:,3], color = :orange, line=(1, [:dash]), label = L"\textrm{ET RV}", marker = ([:x :d], 3, 0.8, stroke(3, :red)))
plot!(dt, VT[:,3], color = :pink, line=(1, [:dash]), label = L"\textrm{VT RV}", marker = ([:+ :d], 3, 0.8, stroke(3, :red)))
hline!([CTHY[3]], color = :brown, line=(1, [:dash]), label = L"\textrm{CT HY}")
hline!([ETHY[3]], color = :lightblue, line=(2, [:dashdot]), label = L"\textrm{ET HY}")
xlabel!(L"\textrm{Sampling interval}")
ylabel!(L"\tilde{\rho}_{\Delta t}^{ij}")
# savefig("Plots/EmpAppCompNEDBTI.svg")

# AGL/NED pair (at index number 5)
plot(dt, CTMM[:,5], legend = :topright, color = :blue, line=(1, [:dash]), label = L"\textrm{CT MM}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)), ylims = (-0.3, 0.3))
plot!(dt, CTRV[:,5], color = :red, line=(1, [:dash]), label = L"\textrm{CT RV}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)))
plot!(dt, ETMM[:,5], color = :green, line=(1, [:dash]), label = L"\textrm{ET MM}", marker = ([:+ :d], 3, 0.8, stroke(3, :blue)))
plot!(dt, ETRV[:,5], color = :orange, line=(1, [:dash]), label = L"\textrm{ET RV}", marker = ([:x :d], 3, 0.8, stroke(3, :red)))
plot!(dt, VT[:,5], color = :pink, line=(1, [:dash]), label = L"\textrm{VT RV}", marker = ([:+ :d], 3, 0.8, stroke(3, :red)))
hline!([CTHY[5]], color = :brown, line=(1, [:dash]), label = L"\textrm{CT HY}")
hline!([ETHY[5]], color = :lightblue, line=(2, [:dashdot]), label = L"\textrm{ET HY}")
xlabel!(L"\textrm{Sampling interval}")
ylabel!(L"\tilde{\rho}_{\Delta t}^{ij}")
# savefig("Plots/EmpAppCompNEDAGL.svg")
