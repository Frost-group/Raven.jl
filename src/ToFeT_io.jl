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

