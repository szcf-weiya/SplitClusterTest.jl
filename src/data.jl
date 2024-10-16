using Copulas

function ar1(ρ::Real, p::Int)
    Σ = zeros(p, p)
    for i = 0:p-1
        for j = 1:p-i
            Σ[j, j+i] = ρ^i
            Σ[j+i, j] = Σ[j, j+i]
        end
    end
    return Σ    
end

function fixcorr_partial(ρ::Real, p::Int; order = 3/4)
    # only keep p^{3/4} correlated terms
    Σ = Matrix(1.0I(p))
    p1 = round(Int, p ^ order)
    idx = sample(1:p, p1, replace=false)
    Σ1 = ones(p1, p1) * ρ + (1 - ρ) * 1.0I
    Σ[idx, idx] .= Σ1
    return Σ
end

# correlation structure: 
#   fixed correlation among S1, 
#   zero correlation among S0, 
#   and fixed correlation between S1 and S0
function fixcorr_s1(ρ::Real, p::Int, p1::Int; corr01 = 1)
    Σ = zeros(p, p)
    Σ[1:p1, 1:p1] .= ρ
    a = min(corr01, sqrt((1 + (p1 - 1) * ρ) / (p1 * (p - p1))) * 0.99) # shrink a little
    Σ[p1+1:p, 1:p1] .= a
    Σ[1:p1, p1+1:p] .= a
    for i = 1:p
        Σ[i, i] = 1
    end
    return Σ
end

function gen_data_normal(n::Int, p::Int, δ = 1.0; prop_imp = 0.1, 
                                corr_structure = "ind", ρ = 0.9,
                                prop_cl1 = 0.5,
                                order = 1,
                                sigma = 0,
                                scale_sigma = 1,
                                two_sided = false, prop_up = 0.5)
    p1 = Int(p * prop_imp)
    if corr_structure == "ar1"
        Σ = ar1(ρ, p)
    elseif corr_structure == "fixcorr_s1"
        Σ = fixcorr_s1(ρ, p, p1, corr01 = 1)
    elseif corr_structure == "fixcorr_s1_ind"
        Σ = fixcorr_s1(ρ, p, p1, corr01 = 0) # Cov(S0, S1) = 0
    elseif corr_structure == "fixcorr"
        if order == 1
            Σ = ones(p, p) * ρ + (1 - ρ) * 1.0I
        else
            Σ = fixcorr_partial(ρ, p, order = order)
        end
    elseif corr_structure == "ind"
        Σ = 1.0I(p)
    else
        @warn "$corr_structure is not supported, use the identity matrix instead"
        Σ = 1.0I(p)
    end
    x = rand(MvNormal(Σ * scale_sigma), n)'
    n1 = Int(n * prop_cl1)
    if two_sided
        for i = 1:p1
            x[1:n1, i] .+= δ * sample([1, -1], Weights([prop_up, 1 - prop_up]))
        end
    else
        x[1:n1, 1:p1] .+= δ
    end
    cl = vcat(ones(Int, n1), ones(Int, n - n1) * 2)
    x = x + randn(n, p) * sigma
    return Matrix(x), cl
end

function gen_data_pois(Λ::AbstractMatrix; ρ = 0.5, block_size = 10)
    n, p = size(Λ)
    Σ = ar1(ρ, block_size)
    cop = GaussianCopula(Σ)
    X = zeros(n, p)
    pp = Int(p / block_size)
    ds = Array{Distribution, 1}(undef, block_size)
    for i = 1:n
        for j = 1:pp
            for k = 1:block_size
                ds[k] = Poisson(Λ[i, (j-1)*block_size+k])
            end
            X[i, (j-1)*block_size+1 : j*block_size] .= rand(SklarDist(cop, ds), 1)
        end
    end
    return X
end

function gen_data_pois(; n = 1000, p = 2000, prop_imp = 0.1, ρ = 0.5, δ = 1.0, block_size = 10, type = "discrete", sigma = 0)
    if type == "discrete"
        L = sample(0:1, n)
    else
        L = randn(n)
        L .-= mean(L)
    end
    p1 = Int(p * prop_imp)
    Fs = vcat(ones(p1)*δ, zeros(p - p1))
    LFT = L * Fs'
    logΛ = LFT .+ log(3) + randn(n, p) * sigma  # intercept
    Λ = exp.(logΛ)
    X = gen_data_pois(Λ, ρ = ρ, block_size = block_size)
    return X
end

