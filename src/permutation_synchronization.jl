function synchronize_permutations{T}(S::Matrix{Matrix{T}})
    m = size(S, 1)
    n = size(S[1], 1)
    Q = Vector{Matrix{Int}}(m)

    W = Symmetric(float(joinblocks(S)))
    U = sqrt(m) * eigvecs(W, ((m-1)*n+1) : (m*n))

    for i in 1:m
        rows = ((i-1)*n+1) : (i*n)
        Pi1 = U[rows, 1:n] * U[1:n, 1:n]'
        Q[i] = eye(Int, n)[:, munkres(-Pi1')]
    end

    return Q
end

function synchronize_permutations{T}(A::Vector{Matrix{T}}, max_depth::Int = 9)
    m = length(A)
    d, n = size(A[1])
    D = register(A[1])
    S = Matrix{Matrix{T}}(m, m)

    for i in 1:m
        S[i, i] = eye(n)

        for j in 1:i-1
            print("\r$i $j")
            flush(STDOUT)

            S[j, i] = register(D, A[i], A[j], max_depth)[2]
            S[i, j] = S[j, i]'
        end
    end

    P = synchronize_permutations(S)

    return [A[i] * P[i] for i in 1:m]
end
