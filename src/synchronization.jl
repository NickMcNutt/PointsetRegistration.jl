import Munkres: munkres

function parallel_pairwise_register{N, T}(n::Int, g::Int, kinds::NTuple{N, Vector{Matrix{T}}}, weights::NTuple{N, T})
    queue = []

    for i in 1:n
        for j in 1:i
            X = ntuple(N) do (k) kinds[k][i] end
            Y = ntuple(N) do (k) kinds[k][j] end
            push!(queue, (weights, X, Y))
        end
    end

    cells = pmap(queue) do pair
        weights, X, Y = pair
        G = grid_coords(2π, g)
        d1, R1, P1 = register(weights, X, Y, G, false)
        d2, R2, P2 = register(weights, X, Y, G, true)
        if d1 < d2
            return d1, R1, P1
        else
            return d2, R2, P2
        end
    end
    
    D = zeros(n, n)
    R = Matrix{Matrix{T}}(n, n)
    P = [Matrix{Matrix{T}}(n, n) for k in 1:N]
    
    for i in 1:n
        for j in 1:i
            d, r, p = shift!(cells)
            D[i, j] = D[j, i] = d
            R[i, j] = r'
            R[j, i] = r

            for k in 1:N
                P[k][i, j] = p[k]'
                P[k][j, i] = p[k]
            end
        end
    end
    
    return D, R, P
end

function project_to_stiefel(C::AbstractMatrix, d::Int)
    n = size(C, 1)
    E = eigfact(Symmetric(C), n-d+1:n)
    return E[:vectors]
end

function stiefel_to_permutations{T}(U::Matrix{T})
    n = size(U, 2)
    m = fld(size(U, 1), n)
    
    P = map(1:m) do i
        P_i1 = U[(i-1)*n + 1:i*n, 1:n] * U[1:n, 1:n]'
        p = munkres(-P_i1')
        eye(T, n)[:, p]
    end
    
    return P
end
