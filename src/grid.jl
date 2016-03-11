function grid{T}(width::T, d::Int)
    x = Array{T}(d, d, d)
    y = Array{T}(d, d, d)
    z = Array{T}(d, d, d)
    
    g = collect(linspace(-width/2, width/2, d))
    
    @inbounds for ix in 1:d, iy in 1:d, iz in 1:d
        x[ix, iy, iz] = g[ix]
        y[ix, iy, iz] = g[iy]
        z[ix, iy, iz] = g[iz]
    end
    
    return x, y, z
end

function grid_coords{T}(width::T, d::Int)
    x, y, z = grid(width, d)
    hcat(vec(x), vec(y), vec(z))'
end
