using LinearAlgebra
using SparseArrays
include("LRA_type.jl")

n = 80
H = Tridiagonal(ones(n-1),-2*ones(n),ones(n-1)) + 0.00001*im*I(n)
H = LRA_mod.LRA(sparse(H),0.0,0.2,0.1)
#display(H[:,:])

δH = 0.01*inv(H)
