const MAX = 10000000

# Lipschitzed versions of f(x) = x^2
f₁(x) = 1 - 1/(1 + x^2)
f₂(x) = 1 - exp(-x^2)

const K₁ = (3/8) * sqrt(3)
const K₂ = sqrt(2) * exp(-1/2)

# Lipschitz (metric) distance between two point sets
function distance_lipschitz{T}(X::Matrix{T}, Y::Matrix{T}, R::Matrix{T}, P::Matrix{Int})
    n = size(Y, 2)

    d = zero(T)
    for i in 1:n
        d += f₁(sumabs2(X[:, i] - (R*Y*P)[:, i]))
    end
    d /= n

    return d
end

# Mean vector length in a point set
function mean_vector_magnitude{T}(Y::Matrix{T})
    n = size(Y, 2)
    
    d = zero(T)
    for i in 1:n
        d += norm(Y[:, i])
    end
    d /= n
    
    return d
end

function eye!{T}(A::Matrix{T})
    for i in eachindex(A)
        A[i] = zero(T)
    end

    for i in 1:min(size(A, 1), size(A, 2))
        A[i, i] = one(T)
    end

    return A
end

immutable RegisterData{T}
    V::Matrix{T}
    W::Matrix{T}
    R::Matrix{T}
    P::Matrix{Int}

    RegisterData{T}(::Type{T}, d::Int, n::Int, MAX::Int = MAX) = new(Matrix{T}(d, MAX), Matrix{T}(d, MAX), Matrix{T}(d, d), Matrix{Int}(n, n))
end

RegisterData{T}(::Type{T}, d::Int, n::Int) = RegisterData{T}(T, d, n)
RegisterData{T}(::Type{T}, d::Int, n::Int, MAX::Int) = RegisterData{T}(T, d, n, MAX)

# Solves the pose and correspondence problem exactly for two point sets
function register{T}(D::RegisterData{T}, X::Matrix{T}, Y::Matrix{T}, max_depth::Int = 6)
    V, W, R, P = D.V, D.W, D.R, D.P

    d, n = size(X)
    L = K₁ * mean_vector_magnitude(Y)
    bw = float(π)
    hbw = bw / 2
    eye!(R)
    eye!(P)
    f_min = distance_lipschitz(X, Y, R, P)

    V[1, 1] = zero(T)
    V[2, 1] = zero(T)
    V[3, 1] = zero(T)
    v = 1

    x_min = zero(T)
    y_min = zero(T)
    z_min = zero(T)

    for depth in 0:max_depth
        w = 0

        for i in 1:v
            sqrt(V[1, i]^2 + V[2, i]^2 + V[3, i]^2) > π && continue

            rotation(R, V[1, i], V[2, i], V[3, i])
            if depth == 0
                eye!(R)
            end

            copy!(P, eye(P)[:, munkres(-X' * R * Y)])
            f = distance_lipschitz(X, Y, R, P)

            if f < f_min
                x_min = V[1, i]
                y_min = V[2, i]
                z_min = V[3, i]
                f_min = f
                #println(f_min)
                #println(vecnorm(X - R*Y*P))
            end

            ϵ = abs(f - f_min) / L
            ϵ > sqrt(d)*hbw && continue

            for sx in (-1, 1), sy in (-1, 1), sz in (-1, 1)
                w += 1
                W[1, w] = V[1, i] + sx*hbw
                W[2, w] = V[2, i] + sy*hbw
                W[3, w] = V[3, i] + sz*hbw
            end
        end

        V, W = W, V
        v = w
        hbw /= 2
    end

    rotation(R, x_min, y_min, z_min)
    copy!(P, eye(P)[:, munkres(-X' * R * Y)])
    nearest_orthogonal(R, X * P' * Y')

    return R, P
end

register{T}(X::Matrix{T}) = RegisterData(T, size(X)...)
register{T}(X::Matrix{T}, MAX::Int) = RegisterData(T, size(X)..., MAX)

# Memory allocating function variants

register{T}(X::Matrix{T}, Y::Matrix{T}) = register(RegisterData(T, size(X)...), X, Y)
register{T}(X::Matrix{T}, Y::Matrix{T}, max_depth::Int) = register(RegisterData(T, size(X)...), X, Y, max_depth)

#println("Hungarian algorithm calls: $t")
#flush(STDOUT)

#println("depth $(depth+1):\t$w / $(8^(depth+1))\n")
#flush(STDOUT)

