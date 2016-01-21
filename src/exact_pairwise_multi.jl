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

immutable RegisterDataMulti{T}
    V::Matrix{T}
    W::Matrix{T}
    R::Matrix{T}
    P::NTuple{4, Matrix{Int}}

    RegisterDataMulti{T}(::Type{T}, d::Int, n::Int, MAX::Int = MAX) = new(Matrix{T}(d, MAX), Matrix{T}(d, MAX), Matrix{T}(d, d),
        (Matrix{Int}(n, n),
         Matrix{Int}(n, n),
         Matrix{Int}(n, n),
         Matrix{Int}(n, n))
        )
end

RegisterDataMulti{T}(::Type{T}, d::Int, n::Int) = RegisterDataMulti{T}(T, d, n)
RegisterDataMulti{T}(::Type{T}, d::Int, n::Int, MAX::Int) = RegisterDataMulti{T}(T, d, n, MAX)

# Solves the pose and correspondence problem exactly for two point sets
function register{T}(D::RegisterDataMulti{T}, X::NTuple{4, Matrix{T}}, Y::NTuple{4, Matrix{T}}, max_depth::Int = 6)
    V, W, R, P = D.V, D.W, D.R, D.P

    d, n = size(X[1])
    L = sum(1:4) do (i) K₁ * mean_vector_magnitude(Y[i]) end

    bw = float(π)
    hbw = bw / 2
    eye!(R)
    foreach(1:4) do (i) eye!(P[i]) end
    f_min = sum(1:4) do (i) distance_lipschitz(X[i], Y[i], R, P[i]) end

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

            foreach(1:4) do (j) copy!(P[j], eye(P[j])[:, munkres(-X[j]' * R * Y[j])]) end
            f = sum(1:4) do (j) distance_lipschitz(X[j], Y[j], R, P[j]) end

            if f < f_min
                x_min = V[1, i]
                y_min = V[2, i]
                z_min = V[3, i]
                f_min = f

                println(f_min)
                s = sum(1:4) do (j) vecnorm(X[j] - R*Y[j]*P[j]) end
                println(s)
                flush(STDOUT)
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

        println("Level $depth")
        flush(STDOUT)
    end

    rotation(R, x_min, y_min, z_min)
    foreach(1:4) do (i) copy!(P[i], eye(P[i])[:, munkres(-X[i]' * R * Y[i])]) end
    XX = hcat(X...)
    YY = hcat(Matrix[Y[i] * P[i] for i in 1:4]...)
    nearest_orthogonal(R, XX * YY')

    return R, P
end

register{T}(X::NTuple{4, Matrix{T}}) = RegisterDataMulti(T, size(X[1])...)
register{T}(X::NTuple{4, Matrix{T}}, MAX::Int) = RegisterDataMulti(T, size(X[1])..., MAX)

# Memory allocating function variants

register{T}(X::NTuple{4, Matrix{T}}, Y::NTuple{4, Matrix{T}}) = register(RegisterDataMulti(T, size(X[1])...), X, Y)
register{T}(X::NTuple{4, Matrix{T}}, Y::NTuple{4, Matrix{T}}, max_depth::Int) = register(RegisterDataMulti(T, size(X[1])...), X, Y, max_depth)
