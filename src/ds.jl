using Clustering
using StatsBase
using RCall
using LinearAlgebra
using Distributions

function calc_τ(ms::AbstractVector, q::Float64 = 0.05, offset::Int = 1)
    ts = sort(abs.(ms))
    ret = [maximum(ms)]
    for t in ts
        curr_fdr = (offset + sum(ms .<= -t)) / max(1, sum(ms .>= t))
        if curr_fdr <= q
            append!(ret, t)
        end
    end
    return minimum(ret)
end

function mirror_stat(d1::AbstractArray{T}, d2::AbstractArray{T}) where T <: AbstractFloat
    return (abs.(d1) + abs.(d2)) .* sign.(d1) .* sign.(d2) * sign(sum(d1 .* d2 ))
end

function cluster_diff(x::AbstractMatrix, cl::AbstractVector{Int64}; method = "sign_pval")
    x1 = x[cl .== 1, :]
    x2 = x[cl .== 2, :]
    p = size(x1, 2)
    # pvals = zeros(p)
    ds = zeros(p)
    for i = 1:p
        if method == "sign_pval"
            pval = rcopy(R"t.test($(x1[:, i]), $(x2[:, i]))$p.value")
            ds[i] = sign(mean(x1[:, i]) - mean(x2[:, i])) * -log(max(eps(), pval))
        elseif method == "tstat"
            t = try 
                rcopy(R"t.test($(x1[:, i]), $(x2[:, i]))$statistic") 
            catch e 
                @warn e
                0 
            end # no enough y (or other error?)
            ds[i] = t
        elseif method == "glm_pois"
            sfit = R"summary(glm($(x[:, i]) ~ as.factor($(cl)), family = 'poisson'))"
            pval = rcopy(R"$sfit$coefficients[2,4]")
            est = rcopy(R"$sfit$coefficients[2,1]")
            # println("pval = $pval, est = $est")
            #ds[i] = sign(est) * -log(pval)
            ds[i] = est
        else
            ds[i] = mean(x1[:, i]) - mean(x2[:, i])
        end
    end
    return ds
end

function calc_sign_pvals(X::AbstractMatrix, t::AbstractVector)
    idx = t .!= -1
    if sum(idx) == 0
        return ones(length(t))
    end
    mt = median(t[idx])
    ii = t[idx] .< mt
    n, p = size(X)
    ssigns = zeros(p)
    pvals = zeros(p)
    for i = 1:p
        mean_diff = mean(X[idx, i][ii]) - mean(X[idx, i][.!ii])
        pval = try rcopy(R"pval_glm"(X[idx, i], t[idx]))
        catch e
            @warn e
            1
        end
        ssigns[i] = sign(mean_diff)
        pvals[i] = -log(max(eps(), pval))
        # spvals[i] = sign(mean_diff) * -log(max(eps(), pval))
    end
    return ssigns .* pvals
end

cl_kmeans(x::AbstractMatrix; kw...) = kmeans(x', 2; kw...).assignments

function ti_pca(xx::AbstractMatrix)
    xx = log.(x .+ 1)
    t = svd(xx .- mean(xx, dims = 1)).U[:, 1]
    return t
end

function ds(x::AbstractMatrix; q = 0.05, signal_measure = "tstat", 
                plt = true, ret_ms = false, 
                type = "discrete", 
                cl_method::Function = cl_kmeans, 
                ti_method::Function = ti_pca, 
                oracle_label = nothing,
                kmeans_whiten = false, Σ = nothing, λ = 1e-6)
    n = size(x, 1)
    n2 = round(Int, n/2)
    idx1 = sample(1:n, n2, replace = false)
    idx2 = setdiff(1:n, idx1)
    x1 = x[idx1, :]
    x2 = x[idx2, :]
    if kmeans_whiten
        p = size(x, 2)
        if isnothing(Σ)
            Σ = est_Σ(x) + λ * 1.0I
            #Σ = cov(x) + λ * 1.0I
        end
        ev = eigen(Σ)
        idx = ev.values .> 0
        xc = x * ev.vectors[idx, :] * diagm(1 ./ sqrt.(ev.values[idx])) * ev.vectors[idx, :]'
        x1c = xc[idx1, :]
        x2c = xc[idx2, :]
    else
        x1c = copy(x1)
        x2c = copy(x2)
    end
    if type == "discrete"
        cl1 = cl_method(x1)
        cl2 = cl_method(x2)
        if !isnothing(oracle_label)
            cl_acc1 = calc_cluster_acc(cl1, oracle_label[idx1])
            cl_acc2 = calc_cluster_acc(cl2, oracle_label[idx2])
            println("cl_acc1 = $cl_acc1, cl_acc2 = $cl_acc2")
        end
        d1 = cluster_diff(x1, cl1; method = signal_measure)
        d2 = cluster_diff(x2, cl2; method = signal_measure)
    else
        t1 = ti_method(x1)
        t2 = ti_method(x2)
        if sum(t1 .* t2) < 0
            t2 = -t2
        end
        d1 = calc_sign_pvals(x1, t1)
        d2 = calc_sign_pvals(x2, t2)
    end
    ms = mirror_stat(d1, d2)
    if ret_ms
        return ms
    end
    τ = calc_τ(ms, q)
    m_sel = findall(ms .> τ)
    if plt
        @info "selected features: " m_sel'
        @info "mirror stat: " ms[m_sel]'
        @info "mirror stat of first 20: " ms[1:20]'
        @info "d1 of first 20: " d1[1:20]'
        @info "d2 of first 20: " d2[1:20]'
        @info "selected features (including negative): " findall(abs.(ms) .> τ)'
        fig = histogram(ms, label = "")
        Plots.vline!(fig, [τ], label = "", lw = 3)
        return fig
    else
        return m_sel
    end
end


function calc_inc_rate(mss::AbstractArray, τs::AbstractArray; sum_in_denom = true, q = 0.05)
    M = length(mss)
    p = length(mss[1])
    inc_rate = zeros(p)
    nempty = 0
    for i = 1:M
        m_sel = findall(mss[i] .> τs[i])
        if sum_in_denom
            inc_rate[m_sel] .+= 1
        else
            inc_rate[m_sel] .+= 1 ./ length(m_sel) # if m_sel is empty, then nothing happens
        end
        if length(m_sel) == 0
            nempty += 1
        end
    end
    if sum_in_denom
        # intuitively, the selected items should appear in most of splits
        # but for sum_in_denom = false, it is not so intuitive to set such a cutoff
        if maximum(inc_rate) <= M * 0.5 # For type I error: if H0, allow at most M*q relevant DE genes
            # An interpretation: without prior knowledge, a feature can be DE or not, only when the select probability larger than 0.5 
            # so only choose features with probability larger than 0.5
            return zeros(p)
        end
    end
    if sum_in_denom
        inc_rate ./= sum(inc_rate)
    else
        inc_rate ./= M
    end
    if nempty >= floor(M * 0.5) # it requires large M to make the proportion accurate
        @info "no relevant features are selected due to too many empty sets"
        return zeros(p)
    end
    return(inc_rate)
end


function sel_inc_rate(inc_rate::AbstractArray{T}; q::Float64 = 0.05, tol = 1e-10, verbose = false) where T <: AbstractFloat
    p = length(inc_rate)
    sort_inc_rate = sort(inc_rate)
    if verbose
        @info "sum of inc_rate: " sum(sort_inc_rate)
        if sum(sort_inc_rate) > 0
            @info "max: " maximum(sort_inc_rate)
            @info "min: " minimum(sort_inc_rate[sort_inc_rate .!= 0])
            # @info sort_inc_rate[sort_inc_rate .> 0]
        end
    end

    curr_rate = 0
    for i = 1:p+1
        if curr_rate > q
            if i <= p
                ntie = sum( abs.(sort_inc_rate[i:end] .- sort_inc_rate[i-1]) .< tol )
            else
                ntie = 0
            end
            ntie0 = sum( abs.(sort_inc_rate[1:i-1] .- sort_inc_rate[i-1]) .< tol )
            if ntie >= ntie0
                return findall(inc_rate .>= sort_inc_rate[i-1])
            else
                return findall(inc_rate .> sort_inc_rate[i-1])
            end
        else
            if i == p+1
                return Int[]
            end
            curr_rate += sort_inc_rate[i]
        end
    end
end

function mds(x::AbstractMatrix; M = 10, q = 0.05, signal_measure = "tstat", plt = false, 
                type = "discrete", 
                qm = 0.05, # dummy to be consistent with mds5
                cl_method = cl_kmeans, 
                ti_method = ti_pca,
                kmeans_whiten = false, kw...)
    if M == 1
        return ds(x; ret_ms = false, plt = false, q = q, signal_measure = signal_measure, type = type, ti_method = ti_method, cl_method = cl_method, kmeans_whiten = kmeans_whiten, kw...)
    end
    mss = [ds(x; ret_ms = true, plt = false, q = q, signal_measure = signal_measure, type = type, ti_method = ti_method, cl_method = cl_method, kmeans_whiten = kmeans_whiten, kw...) for _ in 1:M]
    τs = [calc_τ(_ms, q) for _ms in mss]
    inc_rate = calc_inc_rate(mss, τs, sum_in_denom = true)
    ret = sel_inc_rate(inc_rate)
    return ret
end

