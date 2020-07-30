### A Pluto.jl notebook ###
# v0.11.1

using Markdown
using InteractiveUtils

# ╔═╡ c880bbe4-d252-11ea-299e-f1558843723b
using Revise

# ╔═╡ b5b404ac-d19f-11ea-23ce-ab33f91aa208
begin
	using Pkg
	Pkg.activate("../")
	using Raven
	using LinearAlgebra
	
	using Plots
	gr()
	
	import UnicodePlots
	using IterativeSolvers
end

# ╔═╡ 1b2c0ce4-d1a0-11ea-1c93-a1b0c254c74a
#edges=Raven.readToFeTedges(open("../testdata/10chain.edge","r"))
#edges=Raven.readToFeTedges(open("../testdata/scl.edge","r"))
#edges=Raven.generateLinearChain(20, J=0.01)
edges=Raven.generateNearestNeighbourCube(5)

# ╔═╡ 1e1dd7be-d1a0-11ea-1cc8-c56d0ae6af16
UnicodePlots.spy(edges)

# ╔═╡ 62f8d170-d1a4-11ea-27fd-799898035f6b
siteE=randn(edges.m)./1000

# ╔═╡ 2d27f084-d1a0-11ea-29c9-65aac7f32304
#rates=Raven.ratematrix(edges,Raven.rateAdiabatic)/1e13 # no site energies
rates=Raven.ratematrix(edges,siteE,Raven.rateMarcus)/1e13 # with site energy fluct


# ╔═╡ 307e6efe-d1a0-11ea-059e-4d022439bb2f
# OK, let's look at the SVD (dense) method, as a reference
begin
	if rates.m<1000 # Don't even try this if we have a big problem!
		SVD=Raven.solveMasterDense(rates)
	end
end

# ╔═╡ 351c10fe-d1a0-11ea-3536-6dc718e880b0
# Fiddle till we extract the value we want: the very small singular values, and their right vectors Vt
ρsvd=SVD.Vt[rates.m,:]

# ╔═╡ 434a6b3a-d1a0-11ea-339a-efba5cc9ccb6
sum(ρsvd)

# ╔═╡ 61f5ef98-d1a0-11ea-0700-3182aa28da18
# suspicious number, I wonder if this is...
sqrt(rates.m)

# ╔═╡ ab49bf56-d1a0-11ea-1312-a7547dc76221
r=IterativeSolvers.lobpcg(rates,  true, 1, maxiter=10000)

# ╔═╡ e619cc0c-d1a0-11ea-2173-792ee1914df9
ρ= rates.m* abs.(r.X)/sum(abs.(r.X))

# ╔═╡ 66ce27b8-d1a0-11ea-396a-6b0da7966500
# OK, let's check that we are in a steady state condition, by solving the first order ODE statement
rates.m<1000 && exp(collect(rates))*ρ

# ╔═╡ a0d2c4a2-d1a3-11ea-1fe8-b183f7586657
plot(ρ)

# ╔═╡ 2f365d34-d1aa-11ea-107f-91137d9f43e9
plot(siteE)

# ╔═╡ ea5f86a6-d1a0-11ea-3459-3bf8e7e84e90
# check stationary condition by forming matrix exponential (solution to first order rate equation)
# We cannot do this for real systems, as the matrix exponential is a dense object!
rates.m<1000 && exp(collect(rates))*r.X

# ╔═╡ bf03a888-d252-11ea-0906-3f477c6435bd
Raven.generateNearestNeighbourCube(5)

# ╔═╡ Cell order:
# ╠═c880bbe4-d252-11ea-299e-f1558843723b
# ╠═b5b404ac-d19f-11ea-23ce-ab33f91aa208
# ╠═1b2c0ce4-d1a0-11ea-1c93-a1b0c254c74a
# ╠═1e1dd7be-d1a0-11ea-1cc8-c56d0ae6af16
# ╠═62f8d170-d1a4-11ea-27fd-799898035f6b
# ╠═2d27f084-d1a0-11ea-29c9-65aac7f32304
# ╠═307e6efe-d1a0-11ea-059e-4d022439bb2f
# ╠═351c10fe-d1a0-11ea-3536-6dc718e880b0
# ╠═434a6b3a-d1a0-11ea-339a-efba5cc9ccb6
# ╠═61f5ef98-d1a0-11ea-0700-3182aa28da18
# ╠═66ce27b8-d1a0-11ea-396a-6b0da7966500
# ╠═ab49bf56-d1a0-11ea-1312-a7547dc76221
# ╠═e619cc0c-d1a0-11ea-2173-792ee1914df9
# ╠═a0d2c4a2-d1a3-11ea-1fe8-b183f7586657
# ╠═2f365d34-d1aa-11ea-107f-91137d9f43e9
# ╠═ea5f86a6-d1a0-11ea-3459-3bf8e7e84e90
# ╠═bf03a888-d252-11ea-0906-3f477c6435bd
