# This code generates the 2l+l dimensional real irreducible representations of O(3).
#
# For a function defined on the sphere (f: S² → R) and represented in the basis of the
# real spherical harmonics, these matrix representations rotate/reflect the function
# about an axis in 3D space.
#
# The routines in this file are based upon:
#
#   Rotation Matrices for Real Spherical Harmonics, Direct Determination by Recursion
#   Joseph Ivanic and Klaus Ruedenberg
#   J. Phys. Chem., Vol. 100, No. 15, 1996
#
# The article contains numerous errors, and the "Additions and Corrections"
# — J. Phys. Chem. A, Vol. 102, No. 45, 1998 — also contains at least one error,
# for the matrices generated using Table 2 do not constitute an O(3) group
# homomorphism. A few empirical tests reveal that the term "sqrt(1 - δ(m,-1))" should
# probably be changed to "sqrt(1 + δ(m,-1))" in Table 2 for "V_mm'" and "m < 0".
#
# This file contains an optimized implementation of the (corrected) methods in
# the article. A reference implementation is provided in reference_irreps_O3.jl
#
# These methods are dependently typed on (and compile for) distinct values of l.
# Because of this, the first function call for a new value of l may be slow, but
# subsequent calls should be very fast.
#
# Compilation becomes rather slow beyond about l = 12, but if the goal is to rotate
# high bandwidth-limited functions on the sphere, then there exist better methods that
# do not involve explicitly constructing the entire rotation matrix.

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

function irrep_O3_coeff(l::Int, m::Int, n::Int)
    nu = (l + m) * (l - m)
    nv = (1 + δ(m, 0)) * (l + abs(m) - 1) * (l + abs(m))
    nw = (l - abs(m) - 1) * (l - abs(m))
    d = abs(n) < l ? (l + n) * (l - n) : (2l) * (2l - 1)
    u = sqrt(nu / d)
    v = (1/2) * sqrt(nv / d) * (1 - 2δ(m, 0))
    w = (-1/2) * sqrt(nw / d) * (1 - δ(m, 0))
    
    U = P(l, 0, m, n)
    
    if m > 0
        V = :($(sqrt(1 + δ(m, 1))) * $(P(l, 1, m - 1, n)) - $(1 - δ(m, 1)) * $(P(l, -1, -m + 1, n)))
        W = :($(P(l, 1, m + 1, n)) + $(P(l, -1, -m - 1, n)))
    elseif m < 0
        V = :($(1 - δ(m, -1)) * $(P(l, 1, m + 1, n)) + $(sqrt(1 + δ(m, -1))) * $(P(l, -1, -m - 1, n)))
        W = :($(P(l, 1, m - 1, n)) - $(P(l, -1, -m + 1, n)))
    else
        V = :($(P(l, 1, 1, n)) + $(P(l, -1, -1, n)))
        W = :(0)
    end
    
    ex = :($(u == 0 ? 0 : :($u * ($U))) + $(v == 0 ? 0 : :($v * ($V))) + $(w == 0 ? 0 : :($w * ($W))))
    :(N[$(m + l + 1), $(n + l + 1)] = $(ex))
end

@generated irrep_O3!{L}(N::Matrix, M::Matrix, R::Matrix, ::Type{Val{L}}) = Expr(:block, (irrep_O3_coeff(L, i, j) for i in -L:L, j in -L:L)..., :N)

"""
    irrep_O3!(N, M, R)

Generate an irrep of O(3) of dimension 2l+1 from irrep M of dimension 2l-1 and irrep R of dimension 3, and store it in matrix N.
"""
irrep_O3!{T}(N::Matrix{T}, M::Matrix{T}, R::Matrix{T}) = irrep_O3!(N, M, R, Val{fld(size(N, 1), 2)})

"""
    irrep_O3(M, R)

Generate an irrep of O(3) of dimension 2l+1 from irrep M of dimension 2l-1 and irrep R of dimension 3.
"""
irrep_O3{T}(M::Matrix{T}, R::Matrix{T}) = irrep_O3!(Matrix{T}(size(M, 1) + 2, size(M, 2) + 2), M, R)

"""
    irrep_O3(R, l)

Generate an irrep of O(3) of dimension 2l+1 from irrep R of dimension 3.
"""
irrep_O3{T}(R::Matrix{T}, l::Int) = l == 1 ? R : irrep_O3(irrep_O3(R, l - 1), R)
