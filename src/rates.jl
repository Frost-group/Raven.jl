import LinearAlgebra, IterativeSolvers

# TODO: use Physical constants .jl or whatever for this :^)
const k_B=8.617333262145E-5 # in eV/K
const ħ=6.5821E-16 # eV s ; copied + pasted off Google 

function rateAdiabatic(J)
    2*π/ħ * J^2
end

function rateMarcus(J; T=300, λ=0.5, ΔG=0.0)
    2*π/ħ * J^2 * 1/(sqrt(4*π*λ*k_B*T)) * exp(-(λ+ΔG)^2/(4*π*λ*k_B*T))
end

function ratematrix(edges,ratefn)
    # translate back into the dense COO repr; 
    # saves much time and only generates entries for nonzeros
    As,Bs,Js=SparseArrays.findnz(edges)
    
    rates=ratefn.(Js) # this applies for each element
    
    dim=edges.m # number of sites
    diag=dropdims( - sum(SparseArrays.sparse(As,Bs,rates), dims=2), dims=2 ) 
        # the diagonal of the rate matrix is the negative sum of the column rates
    
    SparseArrays.sparse(vcat(As,1:dim),vcat(Bs,1:dim),vcat(rates,diag) )
end

# None of the below actually works / is tested. Here be dragons.

function solveMasterDense(rates)
    vecs = LinearAlgebra.svd(collect(rates))
 # ... calculate mobility here ...
end

function solveMasterCG(rates)
    IterativeSolvers.cg(rates, SparseArrays.sparsevec(rates*ones(rates.m))) 
end

