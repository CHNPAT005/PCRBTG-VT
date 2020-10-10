## Author: Patrick Chang
# Script file for the Hayashi Yoshida
# Supporting Algorithms are at the start of the script
#  Include:
#           - Kanatani Weight Matrix

#---------------------------------------------------------------------------

using Intervals

#---------------------------------------------------------------------------
### Supporting functions

function Kanatani(t1, t2)
    t1 = filter(!isnan, t1)
    t2 = filter(!isnan, t2)
    L1 = length(t1)
    L2 = length(t2)

    W = zeros(L1-1, L2-1)

    for i in 1:(L1-1)
        I1 = Interval{Open, Closed}(t1[i], t1[i+1])
        for j in 1:(L2-1)
            I2 = Interval{Open, Closed}(t2[j], t2[j+1])
            if !isempty(intersect(I1,I2))
                W[i,j] = 1
            end
        end
    end
    return W
end

#---------------------------------------------------------------------------

function HYcorr(p1, p2, t1, t2)
    t1 = filter(!isnan, t1)
    t2 = filter(!isnan, t2)

    p1 = filter(!isnan, p1)
    p2 = filter(!isnan, p2)

    DiffP1 = diff(log.(p1))
    DiffP2 = diff(log.(p2))

    Den = maximum([t1;t2]) - minimum([t1; t2])

    W = Kanatani(t1, t2)

    Sigma = zeros(Float64, 2, 2)
    Sigma[1, 1] = DiffP1' * DiffP1
    Sigma[2, 2] = DiffP2' * DiffP2
    Sigma[1, 2] = DiffP1' * W * DiffP2
    Sigma[2, 1] = Sigma[1, 2]

    var = diag(Sigma)
    sigma = sqrt.(var)
    rho = Sigma ./ (sigma * sigma')

    return rho, Sigma, var
end

function HYcorrFull(p, t)
    # Get size of matrix
    D = size(p)[2]
    ρ = zeros(D, D)
    # Loop through matrix
    for i in 1:(D-1)
        for j in (i+1):D
            res = HYcorr(p[:,i], p[:,j], t[:,i], t[:,j])
            ρ[i,i] = res[1][1,1]; ρ[j,j] = res[1][2,2]
            ρ[i,j] = res[1][1,2]; ρ[j,i] = res[1][2,1]
        end
    end
    return ρ
end
