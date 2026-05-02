
"""
	(::CommonFace)(f, η1, η2, η3, ξ)

Regularizing coordinate transform for parametrization on the unit triangle.

Common face case.
"""
function (::CommonFace)(f, η1, η2, η3, ξ)

    return (ξ^3) *
           ((η1)^2) *
           (η2) *
           (
               f((1 - ξ, ξ - ξ * η1 + ξ * η1 * η2), (1 - (ξ - ξ * η1 * η2 * η3), ξ - ξ * η1)) +
               f((1 - (ξ - ξ * η1 * η2 * η3), ξ - ξ * η1), (1 - ξ, ξ - ξ * η1 + ξ * η1 * η2)) +
               f((1 - ξ, ξ * η1 * (1 - η2 + η2 * η3)), (1 - (ξ - ξ * η1 * η2), ξ * η1 * (1 - η2))) +
               f((1 - (ξ - ξ * η1 * η2), ξ * η1 * (1 - η2)), (1 - ξ, ξ * η1 * (1 - η2 + η2 * η3))) +
               f((1 - (ξ - ξ * η1 * η2 * η3), ξ * η1 * (1 - η2 * η3)), (1 - ξ, ξ * η1 * (1 - η2))) +
               f((1 - ξ, ξ * η1 * (1 - η2)), (1 - (ξ - ξ * η1 * η2 * η3), ξ * η1 * (1 - η2 * η3)))
           )
end


"""
	(::CommonEdge)(f, η1, η2, η3, ξ)

Regularizing coordinate transform for parametrization on the unit triangle.

Common edge case.
"""
function (::CommonEdge)(f, η1, η2, η3, ξ)

    ξη1 = ξ * η1
    η1η2 = η1 * η2
    η2η3 = η2 * η3
    η1η2η3 = η1η2 * η3
    # ξη2 = ξ * η2
    # ξη3 = ξ * η3

    return (ξ^3) * ((η1)^2) * f((1 - ξ, ξη1 * η3), (1 - ξ * (1 - η1η2), ξη1 * (1 - η2))) +
           (ξ^3) *
           ((η1)^2) *
           (η2) *
           (
               f((1 - ξ, ξη1), (1 - ξ * (1 - η1η2η3), ξη1 * η2 * (1 - η3))) +
               f((1 - ξ * (1 - η1η2), ξη1 * (1 - η2)), (1 - ξ, ξη1 * η2η3)) +
               f((1 - ξ * (1 - η1η2η3), ξη1 * η2 * (1 - η3)), (1 - ξ, ξη1)) +
               f((1 - ξ * (1 - η1η2η3), ξη1 * (1 - η2η3)), (1 - ξ, ξη1 * η2))
           )
end


"""
	(::CommonVertex)(f, η1, η2, η3, ξ)

Regularizing coordinate transform for parametrization on the unit triangle.

Common vertex case.
"""
function (::CommonVertex)(f, η1, η2, η3, ξ)

    ξη1 = ξ * η1
    ξη2 = ξ * η2

    return (ξ^3) * η2 * (f((1 - ξ, ξη1), (1 - ξη2, ξη2 * η3)) + f((1 - ξη2, ξη2 * η3), (1 - ξ, ξη1)))
end


function (::PositiveDistance)(f, η1, η2, η3, ξ)

    u = [η1 * (1 - η2), η2]#jacobian of this transformation is (1-η2)
    v = [η3 * (1 - ξ), ξ]#jacobian of this transformation is (1-ξ)

    return (1 - η2) * (1 - ξ) * (f(u, v))
end



"""
	(::CommonFaceQuad)(f, η1, η2, η3, ξ)

Regularizing coordinate transform for parametrization on the unit square: [0,1]² ↦ Γ.

Common face case.
"""
function (::CommonFaceQuad)(f, η1, η2, η3, ξ)

    ξη1 = ξ * η1 # auxiliary

    mξ = (1 - ξ)# auxiliary
    mξη1 = (1 - ξη1)# auxiliary

    # only 4 different terms occur as argument:
    mξη3 = mξ * η3
    mξη3p = mξη3 + ξ

    mξη1η2  = mξη1 * η2
    mξη1η2p = mξη1η2 + ξη1

    return ξ *
           mξ *
           mξη1 *
           (
               f((mξη3, mξη1η2), (mξη3p, mξη1η2p)) +
               f((mξη1η2, mξη3), (mξη1η2p, mξη3p)) +
               f((mξη3, mξη1η2p), (mξη3p, mξη1η2)) +
               f((mξη1η2, mξη3p), (mξη1η2p, mξη3)) +
               f((mξη3p, mξη1η2), (mξη3, mξη1η2p)) +
               f((mξη1η2p, mξη3), (mξη1η2, mξη3p)) +
               f((mξη3p, mξη1η2p), (mξη3, mξη1η2)) +
               f((mξη1η2p, mξη3p), (mξη1η2, mξη3))
           )
end


"""
	(::CommonEdgeQuad)(f, η1, η2, η3, ξ)

Regularizing coordinate transform for parametrization on the unit square: [0,1]² ↦ Γ.

Common edge case.
"""
function (::CommonEdgeQuad)(f, η1, η2, η3, ξ)

    ξη1 = ξ * η1 # occurs as argument (first two kernels calls)
    ξη2 = ξ * η2 # occurs as argument (last four kernels calls)

    mξ = (1 - ξ) # auxiliary
    mξη1 = (1 - ξη1) # auxiliary
    mξη3 = mξ * η3# occurs as argument (first two kernels calls)
    mξη3p = mξη3 + ξ # occurs as argument (first two kernels calls)

    mξη1η3  = mξη1 * η3# occurs as argument (last four kernels calls)
    mξη1η3p = mξη1η3 + ξη1# occurs as argument (last four kernels calls)

    return (ξ^2) * (
        mξ * (f((mξη3p, ξη2), (mξη3, ξη1)) + f((mξη3, ξη2), (mξη3p, ξη1))) +
        mξη1 * (
            f((mξη1η3p, ξη2), (mξη1η3, ξ)) +
            f((mξη1η3p, ξ), (mξη1η3, ξη2)) +
            f((mξη1η3, ξη2), (mξη1η3p, ξ)) +
            f((mξη1η3, ξ), (mξη1η3p, ξη2))
        )
    )
end


"""
	(::CommonVertexQuad)(f, η1, η2, η3, ξ)

Regularizing coordinate transform for parametrization on the unit square: [0,1]² ↦ Γ.

Common vertex case.
"""
function (::CommonVertexQuad)(f, η1, η2, η3, ξ)

    # only 4 different terms occur as argument (ξ is the fourth):
    ξη1 = ξ * η1
    ξη2 = ξ * η2
    ξη3 = ξ * η3

    return (ξ^3) * (f((ξ, ξη1), (ξη2, ξη3)) + f((ξη1, ξ), (ξη2, ξη3)) + f((ξη1, ξη2), (ξ, ξη3)) + f((ξη1, ξη2), (ξη3, ξ)))
end
