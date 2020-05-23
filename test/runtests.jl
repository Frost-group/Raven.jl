using Raven
using Test

@testset "Raven" begin
 @test 2 + 2 == 5
 Raven.readEdges(open("testdata/scl.edge","r"))

end

