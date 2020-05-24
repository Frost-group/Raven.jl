# ToFeT_io.jl
# Code to read in ToFeT file format .edge and .xyz specifications

import DelimitedFiles
import SparseArrays

function readEdges(f)
    raw=DelimitedFiles.readdlm(f)
    
    # - .+1 renumber indices from 1 rather than zero
    # - double up (vcat) indices for making explicit A->B && implicit B->A
    As=Int.(vcat(raw[:,1],raw[:,2])) .+ 1
    Bs=Int.(vcat(raw[:,2],raw[:,2])) .+ 1

    Js=vcat(raw[:,3],raw[:,3])
    
    SparseArrays.sparse(As,Bs,Js)
end

function ratematrix(edges)
    # translate back into the dense CSC repr; 
    # saves much time and only generates entries for nonzeros
    As,Bs,Js=SparseArrays.findnz(edges)
    
    rates=exp.(-Js)
    
    dim=edges.m # number of sites
    diag=dropdims( - sum(SparseArrays.sparse(As,Bs,rates), dims=2), dims=2 )
    
    SparseArrays.sparse(vcat(As,1:dim),vcat(Bs,1:dim),vcat(rates,diag) )
end

