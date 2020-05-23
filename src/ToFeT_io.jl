# ToFeT_io.jl
# Code to read in ToFeT file format .edge and .xyz specifications

import DelimitedFiles
import SparseArrays

function readEdges(f)
    COO=DelimitedFiles.readdlm(f)
    
    # renumber indices from 1 rather than zero
    As=Int.(vcat(COO[:,1],COO[:,2])) .+ 1
    Bs=Int.(vcat(COO[:,2],COO[:,2])) .+ 1

    Js=vcat(COO[:,3],COO[:,3])
    rates=exp.(-Js)

    ratematrix=SparseArrays.sparse(As,Bs,rates)
end

