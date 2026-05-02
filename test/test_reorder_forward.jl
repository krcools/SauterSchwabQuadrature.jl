using SauterSchwabQuadrature
using StaticArrays
using Test

@testset "Testing reorder_forward" begin
    function repeated_reorder(verts1, verts2, strat)
        for i in 1:1000000
            I, j = SauterSchwabQuadrature.reorder_forward(verts1, verts2, strat)
        end
    end

    ### CommonVertex
    verts1 = [SVector{3,Float64}(rand(3)) for i in 1:3]
    verts2 = [SVector{3,Float64}(rand(3)) for i in 1:3]
    r1 = rand(1:3)
    verts2[2] = verts1[1]
    I1, J1, _, _ = reorder(verts1, verts2, CommonVertex(1))
    I2, J2 = SauterSchwabQuadrature.reorder_forward(verts1, verts2, CommonVertex(1))
    @test prod(I1 .== I2) * prod(J1 .== J2) == 1

    @time repeated_reorder(verts1, verts2, CommonVertex(1))


    ### CommonEdge
    verts2[1] = verts1[3]
    I1, J1, _, _ = reorder(verts1, verts2, CommonEdge(1))
    I2, J2 = SauterSchwabQuadrature.reorder_forward(verts1, verts2, CommonEdge(1))
    @test prod(I1 .== I2) * prod(J1 .== J2) == 1

    @time repeated_reorder(verts1, verts2, CommonEdge(1))

    ### CommonFace
    verts2[3] = verts1[2]
    I1, J1, _, _ = reorder(verts1, verts2, CommonFace(1))
    I2, J2 = SauterSchwabQuadrature.reorder_forward(verts1, verts2, CommonFace(1))
    @test prod(I1 .== I2) * prod(J1 .== J2) == 1

    @time repeated_reorder(verts1, verts2, CommonFace(1))
end
