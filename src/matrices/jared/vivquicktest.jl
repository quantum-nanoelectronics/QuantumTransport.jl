using LinearAlgebra
using SparseArrays
include("LRA_type.jl")
using .LRA_mod

n = 512
H = ComplexF64.(Tridiagonal(-ones(n-1),2*ones(n),-ones(n-1))) #+ 0.00001*im*I(n)
#display(eigvals(H))

#H = rand(n,n)
H₀ = LRA(sparse(H),0.2,0.4,0.1)
show(typeof(H₀))

δH = 0.01*inv(H₀)
Hsum = H₀ + δH

show(typeof(Hsum))