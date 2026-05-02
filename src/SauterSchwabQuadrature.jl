module SauterSchwabQuadrature


# -------- used packages
using LinearAlgebra
using StaticArrays

using TestItems

# -------- exportet parts
# types
export SauterSchwabStrategy
export CommonFace, CommonEdge, CommonVertex, PositiveDistance
export CommonFaceQuad, CommonEdgeQuad, CommonVertexQuad

# functions
export sauterschwab_parameterized, reorder



# -------- included files
include("sauterschwabintegral.jl")
include("pulled_back_integrands.jl")
include("reorder_vertices.jl")

end
