
function rateAdiabatic(J)
    2*pi * J^2
end

function rateMarcus(J; T=300, Lambda=0.5)
    
end

function ratematrix(edges,ratefn)
    # translate back into the dense CSC repr; 
    # saves much time and only generates entries for nonzeros
    As,Bs,Js=SparseArrays.findnz(edges)
    
    rates=ratefn.(Js)
    
    dim=edges.m # number of sites
    diag=dropdims( -sum(SparseArrays.sparse(As,Bs,rates), dims=2), dims=2 )
    
    SparseArrays.sparse(vcat(As,1:dim),vcat(Bs,1:dim),vcat(rates,diag) )
end

