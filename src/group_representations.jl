δ(i::Int, j::Int) = i == j ? 1 : 0

function P(l::Int, i::Int, m::Int, n::Int)
    r = 2
    q = l

    if abs(n) < l
        :(R[$(r + i), $(r + 0)] * M[$(q + m), $(q + n)])
    elseif n == l
        :(R[$(r + i), $(r + 1)] * M[$(q + m), $(q + l - 1)] - R[$(r + i), $(r - 1)] * M[$(q + m), $(q - l + 1)])
    elseif n == -l
        :(R[$(r + i), $(r + 1)] * M[$(q + m), $(q - l + 1)] + R[$(r + i), $(r - 1)] * M[$(q + m), $(q + l - 1)])
    end
end

function real_O3_irrep(l::Int, m::Int, n::Int)
    num_u = (l + m) * (l - m)
    num_v = (1 + δ(m, 0)) * (l + abs(m) - 1) * (l + abs(m))
    num_w = (l - abs(m) - 1) * (l - abs(m))
    den = abs(n) < l ? (l + n) * (l - n) : (2l) * (2l - 1)
    u = sqrt(num_u / den)
    v = (1/2) * sqrt(num_v / den) * (1 - 2δ(m, 0))
    w = (-1/2) * sqrt(num_w / den) * (1 - δ(m, 0))
    
    U = P(l, 0, m, n)
    
    if m == 0
        V = :($(P(l, 1, 1, n)) + $(P(l, -1, -1, n)))
        W = :(0)
    elseif m > 0
        V = :($(sqrt(1 + δ(m, 1))) * $(P(l, 1, m - 1, n)) - $(1 - δ(m, 1)) * $(P(l, -1, -m + 1, n)))
        W = :($(P(l, 1, m + 1, n)) + $(P(l, -1, -m - 1, n)))
    elseif m < 0
        V = :($(1 - δ(m, -1)) * $(P(l, 1, m + 1, n)) + $(sqrt(1 + δ(m, -1))) * $(P(l, -1, -m - 1, n)))
        W = :($(P(l, 1, m - 1, n)) - $(P(l, -1, -m + 1, n)))
    end
    
    ex = :($(u == 0 ? 0 : :($u * ($U))) + $(v == 0 ? 0 : :($v * ($V))) + $(w == 0 ? 0 : :($w * ($W))))
    :(N[$(m + l + 1), $(n + l + 1)] = $(ex))
end

@generated real_O3_irrep!{L}(N::Matrix, M::Matrix, R::Matrix, ::Type{Val{L}}) = Expr(:block, (real_O3_irrep(L, i, j) for i in -L:L, j in -L:L)..., :N)
real_O3_irrep!{T}(N::Matrix{T}, M::Matrix{T}, R::Matrix{T}) = real_O3_irrep!(N, M, R, Val{fld(size(N, 1), 2)})
real_O3_irrep{T}(M::Matrix{T}, R::Matrix{T}) = real_O3_irrep!(Matrix{T}(size(M, 1) + 2, size(M, 1) + 2), M, R)
real_O3_irrep{T}(R::Matrix{T}, l::Int) = l == 1 ? R : real_O3_irrep(real_O3_irrep(R, l - 1), R)
