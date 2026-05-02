using SauterSchwabQuadrature
using Test

using LinearAlgebra
using StaticArrays
using CompScienceMeshes
using JuliaFormatter

# --- testsets
@testset "Testing SauterSchwabQuadrature functionality" begin

    @testset "Triangular " begin

        include("parametric_kernel_generator.jl")
        include("local_space.jl")
        include("numquad.jl")
        include("verificationintegral.jl")

        include("test_reordering.jl")
        include("test_cf_tr.jl")
        include("test_ce_tr.jl")
        include("test_cv_tr.jl")
        include("test_pd_tr.jl")
        include("test_reorder_forward.jl")
    end

    @testset "Quadrilateral " begin

        include("quadrilateral_defs.jl")

        include("test_cf_quad.jl")
        include("test_ce_quad.jl")
        include("test_cv_quad.jl")
    end

    @testset "Test formatting of files" begin
        pkgpath = pkgdir(SauterSchwabQuadrature)   # path of this package including name
        @test format(pkgpath, overwrite=false)     # check whether files are formatted according to the .JuliaFormatter.toml 
    end
end

using TestItemRunner
@run_package_tests
