function neg_XᵀRY!{T}(M::Matrix{T}, X::Matrix{T}, Y::Matrix{T}, R::Matrix{T})
    n = size(X, 2)
    
    @inbounds for i in 1:n, j in 1:n
        s1 = X[1, i] * (R[1, 1] * Y[1, j] + R[1, 2] * Y[2, j] + R[1, 3] * Y[3, j])
        s2 = X[2, i] * (R[2, 1] * Y[1, j] + R[2, 2] * Y[2, j] + R[2, 3] * Y[3, j])
        s3 = X[3, i] * (R[3, 1] * Y[1, j] + R[3, 2] * Y[2, j] + R[3, 3] * Y[3, j])
        
        M[i, j] = -(s1 + s2 + s3)
    end
    
    return nothing
end

function wYPXᵀ!{T}(weight::T, Q::Matrix{T}, X::Matrix{T}, Y::Matrix{T}, p::Vector{Int}) 
    w = weight
    n = size(X, 2)
    
    @inbounds for i in 1:n
        Q[1, 1] += w * Y[1, p[i]] * X[1, i]
        Q[1, 2] += w * Y[1, p[i]] * X[2, i]
        Q[1, 3] += w * Y[1, p[i]] * X[3, i]
        
        Q[2, 1] += w * Y[2, p[i]] * X[1, i]
        Q[2, 2] += w * Y[2, p[i]] * X[2, i]
        Q[2, 3] += w * Y[2, p[i]] * X[3, i]
        
        Q[3, 1] += w * Y[3, p[i]] * X[1, i]
        Q[3, 2] += w * Y[3, p[i]] * X[2, i]
        Q[3, 3] += w * Y[3, p[i]] * X[3, i]
    end
    
    return nothing
end

function cost{T}(M::Matrix{T}, p::Vector{Int})
    n = size(M, 1)
    s = zero(T)
    
    @inbounds for i in 1:n
        s += M[i, p[i]]
    end
    
    return s
end

function cost{T}(weight::T, Q::Matrix{T}, M::Matrix{T}, X::Matrix{T}, Y::Matrix{T}, R::Matrix{T})
    neg_XᵀRY!(M, X, Y, R)
    p = Munkres.munkres(M)
    wYPXᵀ!(weight, Q, X, Y, p)
    
    return weight * cost(M, p)
end

function partial_cost{T}(X::Matrix{T}, Y::Matrix{T})
    n = size(X, 2)
    s = zero(T)
    
    @inbounds for j in 1:n
        s1 = X[1, j] * X[1, j] + Y[1, j] * Y[1, j]
        s2 = X[2, j] * X[2, j] + Y[2, j] * Y[2, j]
        s3 = X[3, j] * X[3, j] + Y[3, j] * Y[3, j]
        
        s += s1 + s2 + s3
    end
    
    return s
end
