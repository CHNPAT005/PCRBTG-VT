## Author: Patrick Chang
# Script file to investigate the Epps effect in event time
# using a Hawkes process to simulate the prices.
# We use the price model by Bacry et al. (2013) and compare
# the MM and HY estimator.

## Preamble

cd("/Users/patrickchang1/PCRBTG-VT")

# Estimators
include("../Functions/Correlation Estimators/Dirichlet/NUFFTcorrDK-FGG.jl")
include("../Functions/Correlation Estimators/HY/HYcorr.jl")

# Hawkes
include("../Functions/Hawkes/Hawkes.jl")
