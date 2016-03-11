function kernel_kmeans(K::Matrix, k::Int, ϵ::Float64 = 1e-5)
    n = size(K, 1)
    C = rand(1:k, n)
    C_next = copy(C)
    
    s_total = Inf
    s_total_next = 0
    
    while true
        for i in 1:n
            k_min = 1
            s_min = Inf

            for kk in 1:k
                s1 = K[i, i]
                
                N_kk = 0
                s2 = 0
                for l in 1:n
                    if C[l] == kk
                        N_kk += 1
                        s2 += K[i, l]
                    end
                end
                s2 = (-2/N_kk) * s2

                s3 = 0
                for l in 1:n, m in 1:n
                    if C[l] == kk && C[m] == kk
                        s3 += K[l, m]
                    end
                end
                s3 = (1/N_kk^2) * s3

                s = s1 + s2 + s3
                if s < s_min
                    s_min = s
                    k_min = kk
                end
            end

            C_next[i] = k_min
            s_total_next += s_min
        end
        
        println(s_total_next)
        
        if abs(s_total_next - s_total) < ϵ
            break
        end
        
        s_total, s_total_next = s_total_next, 0
        C, C_next = C_next, C
    end
    
    return [findin(C_next, i) for i in 1:k]
end

function max_k_cut{T}(D::Matrix{T}, k::Int)
    n = size(D, 1)
    Y = Semidefinite(n, n)
    
    constraints = Constraint[]
    
    for i in 1:n, j in 1:n
        if i == j
            push!(constraints, Y[i, j] == T(1))
        elseif k >= 3
            push!(constraints, Y[i, j] >= -inv(k - 1))
        end
    end

    problem = minimize(trace(D*Y), constraints)
    solve!(problem, MosekSolver())
    
    return Convex.evaluate(Y)
end

function factor_k_cut(Y::Matrix, k::Int)
    #Y = max_k_cut(D, k)
    U, S, V = svd(Y)
    S = diagm(S[1:k])
    W = U[:, 1:k] * sqrt(S)
end
