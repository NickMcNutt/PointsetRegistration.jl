import Blosc: compress, decompress

function parallel_density_grid{T}(C::Matrix{T}, width::T, d::Int, σ::T, ϵ::T = 1e-4, num_workers::Int = nworkers())
    num_points = size(C, 2)
    
    num_points_per_worker = fill(fld(num_points, num_workers), num_workers)
    num_points_per_worker[1:mod(num_points, num_workers)] += 1

    params = map(zip(unshift!(cumsum(num_points_per_worker), 0) + 1, num_points_per_worker)) do r
        C[:, range(r[1], r[2])], width, d, σ, ϵ
    end
    
    ρs_compressed = pmap(params, pids = collect(take(workers(), num_workers))) do p
        Blosc.compress(density_grid(p...), level = 9, shuffle = true, itemsize = sizeof(T))
    end

    ρ = mapreduce(.+, zeros(T, d^3), ρs_compressed) do ρ_compressed
        Blosc.decompress(T, ρ_compressed)
    end

    return reshape(ρ, d, d, d)
end

function density_grid{T}(C::Matrix{T}, width::T, d::Int, σ::T, ϵ::T= 1e-4)
    hw = width / 2
    ρ = zeros(T, d, d, d)

    num_points = size(C, 2)
    
    r_max_sq = max(0, -σ^2 * log(ϵ))
    r_max = sqrt(r_max_sq)
    
    b = floor(Int, d * r_max / width) + 1
    
    g = collect(linspace(-hw, hw, d))
    
    @inbounds for i in 1:num_points
        ix = floor(Int, d * (C[1, i] + hw) / width) + 1
        iy = floor(Int, d * (C[2, i] + hw) / width) + 1
        iz = floor(Int, d * (C[3, i] + hw) / width) + 1
        
        x_min = max(1, ix - b)
        y_min = max(1, iy - b)
        z_min = max(1, iz - b)
        
        x_max = min(d, ix + b)
        y_max = min(d, iy + b)
        z_max = min(d, iz + b)
        
        for jx in x_min:x_max, jy in y_min:y_max, jz in z_min:z_max   
            r_sq = (C[1, i] - g[jx])^2 + (C[2, i] - g[jy])^2 + (C[3, i] - g[jz])^2
            
            if r_sq < r_max_sq
                ρ[jx, jy, jz] += exp(-r_sq / σ^2)
            end
        end
    end
    
    return ρ
end

function density_grid{T}(x::Array{T, 3}, y::Array{T, 3}, z::Array{T, 3}, C::Matrix{T}, σ::T, ϵ::T = 1e-4)
    D = size(x)
    ρ = zeros(T, D)

    num_points = size(C, 2)
    
    R = CartesianRange(size(ρ))
    I_min, I_max = first(R), last(R)
    
    x_spread = x[I_max] - x[I_min]
    y_spread = y[I_max] - y[I_min]
    z_spread = z[I_max] - z[I_min]
    
    r_max_sq = max(0, -σ^2 * log(ϵ))
    r_max = sqrt(r_max_sq)
    
    dx = floor(Int, D[1] * r_max / x_spread) + 1
    dy = floor(Int, D[2] * r_max / y_spread) + 1
    dz = floor(Int, D[3] * r_max / z_spread) + 1

    B = CartesianIndex((dx, dy, dz))
    
    @inbounds for i in 1:num_points
        ix = floor(Int, D[1] * (C[1, i] - x[I_min]) / x_spread) + 1
        iy = floor(Int, D[2] * (C[2, i] - y[I_min]) / y_spread) + 1
        iz = floor(Int, D[3] * (C[3, i] - z[I_min]) / z_spread) + 1
        I = CartesianIndex((ix, iy, iz))
        
        for J in CartesianRange(max(I_min, I - B), min(I_max, I + B))
            r_sq = (C[1, i] - x[J])^2 + (C[2, i] - y[J])^2 + (C[3, i] - z[J])^2
            if r_sq < r_max_sq
                ρ[J] += exp(-r_sq / σ^2)
            end
        end
    end
    
    return ρ
end

function density_grid_stats{T}(width::T, d::Int, σ::T, ϵ::T)
    r_max = sqrt(max(0, -σ^2 * log(ϵ)))
    b = floor(Int, r_max / (width / d)) + 1
    num_blocks = (2b + 1)^3
    
    println("Max radius: $(r_max)")
    println("Blocks per dimension: $b")
    println("Total blocks: $(num_blocks)")
    flush(STDOUT)
end
