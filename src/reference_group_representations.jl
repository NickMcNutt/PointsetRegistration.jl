function R(R₀, l::Int, m::Int, n::Int)
    #println("R: l = $l, m = $m, n = $n")
    if abs(m) > l || abs(n) > l
        0
    elseif l == 1
        R₀(m, n)
    elseif l > 1
        u(l, m, n) * U(R₀, l, m, n) +
        v(l, m, n) * V(R₀, l, m, n) +
        w(l, m, n) * W(R₀, l, m, n)
    end
end

δ(i::Int, j::Int) = i == j ? 1 : 0

function u(l::Int, m::Int, n::Int)
    if abs(n) < l
        sqrt(((l + m) * (l - m)) / ((l + n) * (l - n)))
    elseif abs(n) == l
        sqrt(((l + m) * (l - m)) / ((2l) * (2l - 1)))
    end
end

function v(l::Int, m::Int, n::Int)
    if abs(n) < l
        (1/2) *
        sqrt(((1 + δ(m, 0)) * (l + abs(m) - 1) * (l + abs(m))) / ((l + n) * (l - n))) *
        (1 - 2δ(m, 0))
    elseif abs(n) == l
        (1/2) *
        sqrt(((1 + δ(m, 0)) * (l + abs(m) - 1) * (l + abs(m))) / ((2l) * (2l - 1))) *
        (1 - 2δ(m, 0))
    end
end

function w(l::Int, m::Int, n::Int)
    if abs(n) < l
        (-1/2) *
        sqrt(((l - abs(m) - 1) * (l - abs(m))) / ((l + n) * (l - n))) *
        (1 - δ(m, 0))
    elseif abs(n) == l
        (-1/2) *
        sqrt(((l - abs(m) - 1) * (l - abs(m))) / ((2l) * (2l - 1))) *
        (1 - δ(m, 0))
    end
end

function P(R₀, l::Int, i::Int, m::Int, n::Int)
    #println("P: l = $l, i = $i, m = $m, n = $n")
    if abs(n) < l
        R₀(i, 0) * R(R₀, l - 1, m, n)
    elseif n == l
        R₀(i, 1) * R(R₀, l - 1, m,  l - 1) - R₀(i, -1) * R(R₀, l - 1, m, -l + 1)
    elseif n == -l
        R₀(i, 1) * R(R₀, l - 1, m, -l + 1) + R₀(i, -1) * R(R₀, l - 1, m,  l - 1)
    end
end

function U(R₀, l::Int, m::Int, n::Int)
    #println("U: l = $l, m = $m, n = $n")
    if m == 0
        P(R₀, l, 0, 0, n)
    elseif m > 0
        P(R₀, l, 0, m, n)
    elseif m < 0
        P(R₀, l, 0, m, n)
    end
end

function V(R₀, l::Int, m::Int, n::Int)
    #println("V: l = $l, m = $m, n = $n")
    if m == 0
        P(R₀, l, 1, 1, n) + P(R₀, l, -1, -1, n)
    elseif m > 0
        P(R₀, l, 1, m - 1, n) * sqrt(1 + δ(m, 1)) - P(R₀, l, -1, -m + 1, n) * (1 - δ(m, 1))
    elseif m < 0
        P(R₀, l, 1, m + 1, n) * (1 - δ(m, -1)) + P(R₀, l, -1, -m - 1, n) * sqrt(1 + δ(m, -1))
    end
end

function W(R₀, l::Int, m::Int, n::Int)
    #println("W: l = $l, m = $m, n = $n")
    if m == 0
        0
    elseif m > 0
        P(R₀, l, 1, m + 1, n) + P(R₀, l, -1, -m - 1, n)
    elseif m < 0
        P(R₀, l, 1, m - 1, n) - P(R₀, l, -1, -m + 1, n)
    end
end
