import Munkres: munkres
import Combinatorics: permutations
import MoreMatrices: rotation!

function register{N, T}(weights::NTuple{N, T}, X::NTuple{N, Matrix{T}}, Y::NTuple{N, Matrix{T}}, C::Matrix{T}, invert::Bool = false)
    M = Matrix{T}[Matrix{T}(size(X[i], 2), size(X[i], 2)) for i in 1:N]
    P_min = Matrix{T}[Matrix{T}(size(X[i], 2), size(X[i], 2)) for i in 1:N]
    R = eye(T, 3)
    Q = zeros(T, 3, 3)
    Q_min = similar(Q)
    cost_min = T(Inf)
    
    pc = 0
    @inbounds for i in 1:N
        pc += weights[i] * partial_cost(X[i], Y[i])
    end

    m = size(C, 2)
    @inbounds for i in 1:m
        fill!(Q, 0)
        
        x, y, z = C[1, i], C[2, i], C[3, i]
        
        x^2 + y^2 + z^2 > π^2 && continue
            
        rotation!(R, x, y, z)
        
        if invert
            R[3, 1], R[3, 2], R[3, 3] = -R[3, 1], -R[3, 2], -R[3, 3]
        end
        
        # Q += (wYPXᵀ)₁ + (wYPXᵀ)₂ + ...
        for j in 1:N
            weights[j] != 0 && cost(weights[j], Q, M[j], X[j], Y[j], R)
        end
        
        # Refine
        nearest_orthogonal!(R, Q')

        fill!(Q, 0)
        # Q += (wYPXᵀ)₁ + (wYPXᵀ)₂ + ...
        for j in 1:N
            weights[j] != 0 && cost(weights[j], Q, M[j], X[j], Y[j], R)
        end

        # End refine

        c = pc - 2 * sum(svdvals(Q))
        
        if c < cost_min
            cost_min = c

            copy!(Q_min, nearest_orthogonal(Q))

            for j in 1:N
                p = munkres(M[j])
                P_min[j] = eye(M[j])[:, p]
            end
        end
    end
    
    return sqrt(abs(cost_min)), Q_min, P_min
end

function exact_register{T}(X::Matrix{T}, Y::Matrix{T})
    d, n = size(X)
    s_max = -Inf
    p_min = Vector{Int}(n)
    for p in permutations(1:n)
        s = sum(svdvals(Y[:, p] * X'))
        if s > s_max
            s_max = s
            p_min = p
        end
    end
    R_min = nearest_orthogonal(Y[:, p_min] * X')'
    P_min = eye(n)[:, p_min]

    return vecnorm(X - R_min * Y * P_min), R_min, P_min
end

function gaussian_register{T}(σ::T, X::Matrix{T}, Y::Matrix{T}, G::Matrix{T})
    num_samples = size(G, 2)
    n = size(X, 2)
    R = Matrix{T}(3, 3)
    R_min = Matrix{T}(3, 3)
    
    s_max = 0
    for k in 1:num_samples
        x, y, z = G[1, k], G[2, k], G[3, k]
        
        x^2 + y^2 + z^2 > π^2 && continue
        
        rotation!(R, x, y, z)
        
        s = 0
        for i in 1:n, j in 1:n
            s += exp(vecnorm(X[:, i] - R * Y[:, j])^2 / (-2σ^2))
        end
        
        if s > s_max
            s_max = s
            copy!(R_min, R)
        end
    end
    
    return R_min
end

function register_to_reference{N, T}(g::Int, reference::NTuple{N, Matrix{T}}, kinds::NTuple{N, Vector{Matrix{T}}}, weights::NTuple{N, T})
    n = length(kinds[1])
    queue = []
    
    for i in 1:n
        X = ntuple(N) do (k) reference[k] end
        Y = ntuple(N) do (k) kinds[k][i] end
        push!(queue, (weights, X, Y))
    end
        
    cells = pmap(queue) do pair
        weights, X, Y = pair
        G = grid_coords(1π, g)   
        d1, R1, P1 = register(weights, X, Y, G, false)
        d2, R2, P2 = register(weights, X, Y, G, true)
        if d1 < d2
            return d1, R1, P1
        else
            return d2, R2, P2
        end
    end
    
    R = Vector{Matrix{T}}(n)
    P = ntuple(N) do (k) Vector{Matrix{T}}(n) end
    
    for i in 1:n
        d, r, p = shift!(cells)
        R[i] = r
        
        for k in 1:N
            P[k][i] = p[k]'
        end
    end
    
    return R, P 
end
