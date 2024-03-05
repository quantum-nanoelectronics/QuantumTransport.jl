using LinearAlgebra
using SparseArrays
include("LRA_type.jl")

n = 512
H = ComplexF64.(Tridiagonal(-ones(n-1),2*ones(n),-ones(n-1))) #+ 0.00001*im*I(n)
#display(eigvals(H))

#H = rand(n,n)
H = LRA_mod.LRA(sparse(H),0.4,0.8,0.1)
show(H)

#δH = 0.01*inv(H)
