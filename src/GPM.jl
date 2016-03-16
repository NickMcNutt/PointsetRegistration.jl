import Munkres: munkres

function GPM_orthogonal{T}(C::Matrix{T}, x₁::Matrix{T}, num_iters::Int, r::Int = 3, c::Int = 3)
    α = max(0, -eigmin(C))
    C̃ = C + α*I
    
    x = Vector{Matrix{T}}(num_iters)
    x[1] = x₁
    
    for i in 2:num_iters
        D = C̃ * x[i-1]
        O = split_blocks(D, r, c)
        map!(nearest_orthogonal, O)
        x[i] = join_blocks(O)
    end
    
    return x
end

function GPM_permutation{T}(C::Matrix{T}, x₁::Matrix{T}, num_iters::Int, d::Int)
    α = max(0, -eigmin(C))
    C̃ = C + α*I
    
    x = Vector{Matrix{T}}(num_iters)
    x[1] = x₁
    
    for i in 2:num_iters
        D = C̃ * x[i-1]
        O = split_blocks(D, d, d)
        map!(O) do Q
            eye(Int, d)[:, munkres(-Q')]
        end
        x[i] = join_blocks(O)
    end
    
    return x
end
