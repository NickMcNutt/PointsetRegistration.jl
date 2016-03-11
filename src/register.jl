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
            
        MoreMatrices.rotation!(R, x, y, z)
        
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
                p = Munkres.munkres(M[j])
                P_min[j] = eye(M[j])[:, p]
            end
        end
    end
    
    return sqrt(abs(cost_min)), Q_min, P_min
end

function iterative_register{N, T}(w::NTuple{N, T}, S::Vector{Matrix{T}}, A::Vector{Matrix{T}}, G::Matrix{T}, invert::Bool = false)
    s_max = 0
    R_max = Matrix{T}(3, 3)
    p_max = Vector{Vector{Int}}(N)
    
    R = eye(T, 3)
    p = Vector{Vector{Int}}(N)
    Q = Matrix{T}(3, 3)
    
    m = size(G, 2)
    for i in 1:m
        x, y, z = G[1, i], G[2, i], G[3, i]
        
        x^2 + y^2 + z^2 > π^2 && continue
            
        MoreMatrices.rotation!(R, x, y, z)
        
        if invert
            R[3, 1], R[3, 2], R[3, 3] = -R[3, 1], -R[3, 2], -R[3, 3]
        end   
        
        fill!(Q, 0)
        
        for j in 1:N
            p[j] = munkres(-(S[j] * R * A[j])')
            Q += w[j] * A[j][:, p[j]] * S[j]
        end
        
        s = sum(svdvals(Q))
        
        if s > s_max
            s_max = s
            copy!(R_max, nearest_orthogonal(Q'))
            
            for j in 1:N
                p_max[j] = p[j]
            end
        end
    end
    
    return s_max, R_max, p_max
end

function iterative_register{N, T}(weights::NTuple{N, T}, X::NTuple{N, Vector{Matrix{T}}})
    m = length(X[1])
    G = grid_coords(1π, 16)
    A = Vector{Matrix{T}}[deepcopy(Y) for Y in X]
    B = Vector{Matrix{T}}(N)
    S = Vector{Matrix{T}}(N)
    
    @inbounds for t in 1:8
        for k in 1:m
            s_max = 0
            R_max = Matrix{T}(3, 3)
            p_max = Vector{Vector{Int}}(N)

            for j in 1:N
                S[j] = sum(A[j])' - A[j][k]'
                S[j] -= A[j][k]'
                B[j] = A[j][k]
            end     
            
            s_max, R_max, p_max = iterative_register(weights, S, B, G, false)
            s_max2, R_max2, p_max2 = iterative_register(weights, S, B, G, true)
            
            if s_max2 > s_max
                s_max = s_max2
                R_max = R_max2
                p_max = p_max2
            end
            
            for j in 1:N
                A[j][k] = R_max * B[j][:, p_max[j]]
            end
            
            s = 0
            for j in 1:N
                s += weights[j] * vecnorm(mean(A[j]))^2
            end
                
            print("\r$t\t$k: $s")
            flush(STDOUT)
        end
        
        println()
    end
    
    return map(mean, A), A
end

function register_exact{T}(X::Matrix{T}, Y::Matrix{T})
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
