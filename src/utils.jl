"""
    calc_acc(pred::AbstractVector, truth::AbstractVector)

Calculate the accuracy (FDR, F1 score, Precision, Power) for the predicted selection `pred` given the `truth`.
"""
function calc_acc(pred::AbstractVector{T}, truth::AbstractVector{T}) where T <: Real
    if length(truth) == 0 # NULL
        if length(pred) == 0
            return 0, 1, 1, 1
        else
            return 1, 0, 0, 0            
        end
    end
    if length(pred) == 0 # truth not empty
        return 0, 0, 0, 0
    end
    tp = length(intersect(pred, truth))
    # if tp == 0
    #     return 1, 0, 0, 0
    # end
    recall = tp / length(truth)
    precision = tp / length(pred)
    fdr = 1 - precision
    f1 = 2 * precision * recall / (precision + recall)
    if isnan(f1)
        f1 = 0
    end
    return fdr, f1, precision, recall
end

"""
    first_two_PCs(x::AbstractMatrix)

Calculate the first two principal components of data matrix `x`.
"""
function first_two_PCs(x::AbstractMatrix)
    xc = x .- mean(x, dims = 1)
    U, S, V = svd(xc)
    return U[:, 1] * S[1], U[:, 2] * S[2]
end

function est_Σ(x::AbstractMatrix)
    n, p = size(x)
    cl = rcopy(R"kmeans($x, 2, nstart = 10)$cluster") # the accuracy is low
    Σ1 = cov(x[cl .== 1, :])
    Σ2 = cov(x[cl .== 2, :])
    Σ = (Σ1 * sum(cl .== 1) + Σ2 * sum(cl .== 2)) / (n - 2)
    return Σ
end

