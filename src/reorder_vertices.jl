function reorder(t, s, strat::CommonVertex)

    T = eltype(t[1])
    tol = 1e3 * eps(T)
    # tol = 1e5 * eps(T)
    # tol = sqrt(eps(T))

    # Find the permutation P of t and s that make
    # Pt = [P, A1, A2]
    # Ps = [P, B1, B2]
    I = zeros(Int, 1)
    J = zeros(Int, 1)
    e = 1
    for i in 1:3
        v = t[i]
        for j in 1:3
            w = s[j]
            if norm(w - v) < tol
                I[e] = i
                J[e] = j
                e += 1
                break
            end
        end
        e == 2 && break
    end

    append!(I, setdiff([1, 2, 3], I))
    append!(J, setdiff([1, 2, 3], J))

    # # inverse permutations
    # K = indexin([1,2,3], I)
    # L = indexin([1,2,3], J)

    K = zeros(Int, 3)
    for i in 1:3
        for j in 1:3
            if I[j] == i
                K[i] = j
                break
            end
        end
    end

    L = zeros(Int, 3)
    for i in 1:3
        for j in 1:3
            if J[j] == i
                L[i] = j
                break
            end
        end
    end

    return I, J, K, L
end


function reorder(t, s, strat::CommonEdge)

    T = eltype(t[1])
    tol = 1e3 * eps(T)
    # tol = 1e5 * eps(T)
    # tol = sqrt(eps(T))

    I = zeros(Int, 3)
    J = zeros(Int, 3)
    e = 1
    for i in 1:3
        v = t[i]
        for j in 1:3
            w = s[j]
            if norm(w - v) < tol
                I[e] = i
                J[e] = j
                e += 1
                break
            end
        end
    end
    I[3] = setdiff([1, 2, 3], I[1:2])[1]
    J[3] = setdiff([1, 2, 3], J[1:2])[1]

    I = circshift(I, -1)
    J = circshift(J, -1)

    # # inverse permutations
    # K = indexin([1,2,3], I)
    # L = indexin([1,2,3], J)

    K = zeros(Int, 3)
    for i in 1:3
        for j in 1:3
            if I[j] == i
                K[i] = j
                break
            end
        end
    end

    L = zeros(Int, 3)
    for i in 1:3
        for j in 1:3
            if J[j] == i
                L[i] = j
                break
            end
        end
    end

    return I, J, K, L
end


function reorder(t, s, strat::CommonFace)

    T = eltype(t[1])
    tol = 1e3 * eps(T)
    # tol = 1e5 * eps(T)
    # tol = sqrt(eps(T))


    I = [1, 2, 3]
    J = [-1, -1, -1]
    numhits = 0
    for (i, v) in pairs(t)
        for (j, w) in pairs(s)
            if norm(w - v) < tol
                J[i] = j
                numhits += 1
            end
        end
    end

    @assert numhits == 3
    @assert all(J .!= -1)

    K = zeros(Int, 3)
    for i in 1:3
        for j in 1:3
            if I[j] == i
                K[i] = j
                break
            end
        end
    end

    L = zeros(Int, 3)
    for i in 1:3
        for j in 1:3
            if J[j] == i
                L[i] = j
                break
            end
        end
    end

    return I, J, K, L
end


# Find the permutation P of t and s that make
# Pt = [P, A1, A2, A3]
# Ps = [P, B1, B2, B3]
function reorder(t, s, strat::CommonVertexQuad)
    T = eltype(eltype(t))
    T = eltype(t[1])
    tol = 1e3 * eps(T)

    I = zeros(Int, 1)
    J = zeros(Int, 1)
    e = 1
    for i in 1:4
        v = t[i]
        for j in 1:4
            w = s[j]
            if norm(w - v) < tol
                I[e] = i
                J[e] = j
                e += 1
                break
            end
        end
        e == 2 && break
    end

    I = circshift([1, 2, 3, 4], 1 - I[1])
    J = circshift([1, 2, 3, 4], 1 - J[1])

    return I, J, nothing, nothing
end


@testitem "reorder CommonVertexQuad" begin
    t = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]

    s = [[0.0, -1.0, 0.0], [-1.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

    quadPW = SauterSchwabQuadrature._legendre(10, 0, 1)
    strat = SauterSchwabQuadrature.CommonVertexQuad(quadPW)
    I, J, K, L = SauterSchwabQuadrature.reorder(t, s, strat)

    @show I
    @show J

    @test I[1] == 1
    @test J[1] == 4
end


function reorder(t, s, strat::CommonEdgeQuad)
    T = eltype(t[1])
    tol = 1e3 * eps(T)

    I = zeros(Int, 2)
    J = zeros(Int, 2)
    e = 1
    for i in 1:4
        v = t[i]
        for j in 1:4
            w = s[j]
            if norm(w - v) < tol
                I[e] = i
                J[e] = j
                e += 1
                break
            end
        end
        e == 3 && break
    end

    if mod1(I[1] + 1, 4) == I[2]
        I = circshift([1, 2, 3, 4], 1 - I[1])
    else
        I = circshift([4, 3, 2, 1], 1 - (5 - I[1]))
    end

    if mod1(J[1] + 1, 4) == J[2]
        J = circshift([1, 2, 3, 4], 1 - J[1])
    else
        J = circshift([4, 3, 2, 1], 1 - (5 - J[1]))
    end

    # append!(I, setdiff([1, 2, 3, 4], I))
    # append!(J, setdiff([1, 2, 3, 4], J))

    return I, J, nothing, nothing
end

@testitem "reorder CommonEdgeQuad" begin
    t = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]

    s = [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

    quadPW = SauterSchwabQuadrature._legendre(10, 0, 1)
    strat = SauterSchwabQuadrature.CommonEdgeQuad(quadPW)
    I, J, _, _ = SauterSchwabQuadrature.reorder(t, s, strat)

    @show I
    @show J

    @test t[I[1]] ≈ s[J[1]]
    @test t[I[2]] ≈ s[J[2]]
end


function reorder(t, s, strat::CommonFaceQuad)
    T = eltype(eltype(t))
    tol = 1e3 * eps(T)

    I = [1, 2, 3, 4]
    J = [-1, -1, -1, -1]
    numhits = 0
    for (i, v) in pairs(t)
        for (j, w) in pairs(s)
            if norm(w - v) < tol
                J[i] = j
                numhits += 1
            end
        end
    end

    @assert numhits == 4
    @assert all(J .!= -1)

    return I, J, nothing, nothing
end


@testitem "reorder CommonFaceQuad" begin
    t = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]

    s = [[1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]

    quadPW = SauterSchwabQuadrature._legendre(10, 0, 1)
    strat = SauterSchwabQuadrature.CommonFaceQuad(quadPW)
    I, J, _, _ = SauterSchwabQuadrature.reorder(t, s, strat)

    @show I
    @show J

    @test J == [3, 4, 1, 2]
end



#implements allocation free reorder routines without accesing the inverse permutations K and L

function genI(n)
    if n == 1
        return (1, 2, 3)
    elseif n == 2
        return (2, 1, 3)
    elseif n == 3
        return (3, 1, 2)
    end
end

function reorder_forward(t, s, strat::CommonVertex)

    T = eltype(t[1])
    tol = 1e3 * eps(T)

    I = 1
    J = 1
    e = 1
    for i in 1:3
        v = t[i]
        for j in 1:3
            w = s[j]
            if norm(w - v) < tol
                I = i
                J = j
                e += 1
                break
            end
        end
        e == 2 && break
    end


    return genI(I), genI(J)
end

function genI(i1, i2)
    if i1 == 1
        if i2 == 2
            return (2, 3, 1)
        elseif i2 == 3
            return (3, 2, 1)
        end
    elseif i1 == 2
        if i2 == 3
            return (3, 1, 2)
        elseif i2 == 1
            return (1, 3, 2)
        end
    elseif i1 == 3
        if i2 == 1
            return (1, 2, 3)
        elseif i2 == 2
            return (2, 1, 3)
        end
    end
    return (1,2,3)
end

function reorder_forward(t, s, strat::CommonEdge)

    T = eltype(t[1])
    tol = 1e3 * eps(T)
    i1, i2 = 1, 2
    j1, j2 = 1, 2
    e = 1
    for i in 1:3
        v = t[i]
        for j in 1:3
            w = s[j]
            if norm(w - v) < tol
                if e == 1
                    i1 = i
                    j1 = j
                elseif e == 2
                    i2 = i
                    j2 = j
                end
                e += 1
                break
            end
            # e == 3 && break
        end
        e == 3 && break
    end
    return genI(i1, i2), genI(j1, j2)
end
function assign(i1, i2, i3, j1, j2, j3)
    if (i1, i2, i3) == (1, 2, 3)
        return (j1, j2, j3)
    elseif (i1, i2, i3) == (2, 3, 1)
        return (j2, j3, j1)
    elseif (i1, i2, i3) == (3, 1, 2)
        return (j3, j1, j2)
    elseif (i1, i2, i3) == (1, 3, 2)
        return (j1, j3, j2)
    elseif (i1, i2, i3) == (2, 1, 3)
        return (j2, j1, j3)
    elseif (i1, i2, i3) == (3, 2, 1)
        return (j3, j2, j1)
    end
end

function reorder_forward(t, s, strat::CommonFace)

    T = eltype(t[1])
    tol = 1e3 * eps(T)
    # tol = 1e5 * eps(T)
    # tol = sqrt(eps(T))


    I = (1, 2, 3)
    j1, j2, j3 = 1, 2, 3
    i1, i2, i3 = 1, 2, 3
    e = 1
    # numhits = 0
    for (i, v) in pairs(t)
        for (j, w) in pairs(s)
            if norm(w - v) < tol
                if e == 1
                    j1 = j
                    i1 = i
                elseif e == 2
                    j2 = j
                    i2 = i
                elseif e == 3
                    j3 = j
                    i3 = i
                end
                e += 1
            end
            e == 4 && break
        end
        e == 4 && break
    end

    # @assert numhits == 3
    # @assert all(J .!= -1)

    return I, assign(i1, i2, i3, j1, j2, j3)
end
reorder_forward(a, b, c) = reorder(a, b, c)[1:2]

