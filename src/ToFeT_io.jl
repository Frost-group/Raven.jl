# ToFeT_io.jl
# Code to read in ToFeT file format .edge and .xyz specifications

import DelimitedFiles
import SparseArrays

function readToFeTedges(f)
    raw=DelimitedFiles.readdlm(f)

    # - .+1 renumber indices from 1 rather than zero
    # - double up (vcat) indices for making explicit A->B && implicit B->A
    As=Int.(vcat(raw[:,1],raw[:,2])) .+ 1
    Bs=Int.(vcat(raw[:,2],raw[:,1])) .+ 1

    Js=vcat(raw[:,3],raw[:,3])
    
    SparseArrays.sparse(As,Bs,Js)
end

function generateLinearChain(n; J=0.1)
    # bit ugly, surely a nicer way to write this?
    As=[ i for i=1:n-1 ]
    Bs=[ i+1 for i=1:n-1 ]
    Js=[J for i=1:n-1]

    SparseArrays.sparse(vcat(As,Bs),vcat(Bs,As),vcat(Js,Js))
end

function generateNearestNeighbourCube(n; J=0.1)
    N=n^3 # number of sites in cube
    As=Int[]
    Bs=Int[]
    Js=Float64[]

    nn=[+1 -1 +n -n +n^2 -n^2]

    for i in 1:N # iterate over all elements
        for neighbour in nn
            append!(As,i)
            append!(Bs,(i+neighbour-1+N)%N + 1) # modulo arithmatic to enforce PBCs
            append!(Js,J)
        end
    end

    SparseArrays.sparse(As,Bs,Js)
end

function readToFeTxyz(f)
     raw=DelimitedFiles.readdlm(f)
     XYZ=raw[:,1:3]
 end

