using MoreMatrices

trim{T}(coords::Vector{Matrix{T}}, n::Int, m::Int) = [coords[i][:, 1:n] for i in 1:m]

function random_perturb!{T}(A::Vector{Matrix{T}}, σ::T)
    m = length(A)
    d, n = size(A[1])
    N = Vector{Matrix{T}}(m)

    for i in 1:m
        N[i] = σ * randn(d, n)
        copy!(A[i], A[i] + N[i])
    end

    return -N
end

function random_permute!{T}(A::Vector{Matrix{T}})
    m = length(A)
    d, n = size(A[1])
    P = Vector{Matrix{Int}}(m)

    for i in 1:m
        P[i] = eye(n)[:, randperm(n)]
        copy!(A[i], A[i] * P[i]')
    end

    return P
end

function random_rotate!{T}(A::Vector{Matrix{T}})
    m = length(A)
    d, n = size(A[1])
    R = Vector{Matrix{T}}(m)

    for i in 1:m
        R[i] = random_orthogonal(3)
        copy!(A[i], R[i]' * A[i])
    end

    return R
end



function pointsets_transform{T}(A::Vector{Matrix{T}}, R::Vector{Matrix{T}}, P::Vector{Matrix{T}})
    d, n = size(A[1])
    m = min(length(R), length(A), length(P))
    A_hat = [R[i] * A[i] * P[i] for i in 1:m]
end

function pointsets_transform{T}(A::Vector{Matrix{T}}, P::Vector{Matrix{T}})
    d, n = size(A[1])
    m = min(length(R), length(A), length(P))
    A_hat = [A[i] * P[i] for i in 1:m]
end 

function pointsets_transform{T}(A::Vector{Matrix{T}}, R::Vector{Matrix{T}})
    d, n = size(A[1])
    m = min(length(R), length(A), length(P))
    A_hat = [R[i] * A[i] for i in 1:m]
end



function pointsets_error{T}(A::Vector{Matrix{T}})
    d, n = size(A[1])
    m = length(A)
    A_hat = [A[i] for i in 1:m]
    X_hat = mean(A_hat)
    return vecnorm([vecnorm(X_hat - A_hat[i]) for i in 1:m]) / sqrt(d*n*m)
end

function pointsets_error{T}(A::Vector{Matrix{T}}, R::Vector{Matrix{T}})
    d, n = size(A[1])
    m = length(A)
    A_hat = [R[i] * A[i] for i in 1:m]
    X_hat = mean(A_hat)
    return vecnorm([vecnorm(X_hat - A_hat[i]) for i in 1:m]) / sqrt(d*n*m)
end

function pointsets_error{T}(A::Vector{Matrix{T}}, R::Vector{Matrix{T}}, P::Vector{Matrix{Int}})
    d, n = size(A[1])
    m = length(A)
    A_hat = [R[i] * A[i] * P[i] for i in 1:m]
    X_hat = mean(A_hat)
    return vecnorm([vecnorm(X_hat - A_hat[i]) for i in 1:m]) / sqrt(d*n*m)
end

function pointsets_error{T}(A::Vector{Matrix{T}}, R::Vector{Matrix{T}}, P::Vector{Matrix{T}})
    d, n = size(A[1])
    m = length(A)
    A_hat = [R[i] * A[i] * P[i] for i in 1:m]
    X_hat = mean(A_hat)
    return vecnorm([vecnorm(X_hat - A_hat[i]) for i in 1:m]) / sqrt(d*n*m)
end
