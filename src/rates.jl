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

function bareratematrix(edges,ratefn)
    # translate back into the dense COO repr; 
    # saves much time and only generates entries for nonzeros
    As,Bs,Js=SparseArrays.findnz(edges)
    
    rates=ratefn.(Js) # this applies for each element
    
    SparseArrays.sparse(As, Bs, rates )
end

"Single iteration of Yu's power iterative method to solve master equation. Needs a bare ratematrix. (i.e. zeros / empty on the diagonal).

See IV.C. (equation 20):
Yu, Z. G., Smith, D. L., Saxena, A., Martin, R. L., Bishop, A. R. (2001). 
Molecular geometry fluctuations and field-dependent mobility in conjugated polymers. 
Physical Review B, 63(8), 085202. https://doi.org/10.1103/PhysRevB.63.085202
"
function YuPRBpowm!(ω, P)
    Σ=sum # just for notation niceness
    sP=SparseArrays.sparse(P) # sparse representation
	⋅=SparseArrays.dot # dot product for sum reduce
	for i in eachindex(P)
        P[i] = ω[i,:]⋅sP / Σ(ω[:,i]) / (1 - (ω[:,i].-ω[i,:])⋅sP / Σ(ω[:,i]) ) 
    end
    P
end

"Primitive loop to run YuPRBpowm to convergence."
function YuPRBpowm_converge!(rates,P; eps=0.0001)
    c=0
    while true
        c=c+1
        Pold=copy(P)
        Raven.YuPRBpowm!(rates,P)
   		@show c 
        if sum(abs.(P-Pold))<eps
            break
        end
    end
    return c
end

# None of the below actually works / is tested. Here be dragons.

function solveMasterDense(rates)
    vecs = LinearAlgebra.svd(collect(rates))
 # ... calculate mobility here ...
end

function solveMasterCG(rates)
    IterativeSolvers.cg(rates, SparseArrays.sparsevec(rates*ones(rates.m))) 
end

