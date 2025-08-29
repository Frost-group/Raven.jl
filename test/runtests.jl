using Raven
using Test
using SparseArrays

@testset "Raven tests" begin
    # Test testing 
    @test 2 + 2 == 4

    # Test reading edge file
    edges = Raven.readToFeTedges(open("../testdata/10chain.edge","r"))
    @test typeof(edges) == SparseArrays.SparseMatrixCSC{Float64, Int64}
    @test size(edges) == (10, 10)  # 10 sites in the chain
    
    # Test rate functions positive definite
    @test Raven.rateAdiabatic(0.1) > 0.0
    @test Raven.rateMarcus(0.1) > 0.0
    
    # Test rate matrix generation
    rate_matrix = Raven.ratematrix(edges, Raven.rateAdiabatic)
    @test typeof(rate_matrix) == SparseArrays.SparseMatrixCSC{Float64, Int64}
    @test size(rate_matrix) == (10, 10)
end

