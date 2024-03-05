using LinearAlgebra
using SparseArrays
include("LowRankMatrices.jl")
using .LowRankMatrices
#import .LowRankMatrices

n = 8
H = ComplexF64.(Tridiagonal(-ones(n-1),1*ones(n),-ones(n-1))) + 0.00001*im*I(n)
#display(eigvals(H))

#H = rand(n,n)
H₀ = LRA(sparse(H), -0.2, 0.2, initialrank=4, rankstep=2)

show(H₀)

show(typeof(H₀))

δH = 1E-3*inv(H₀)
show(δH)
display(real.(Array(H₀)))
display(Array(H))
#Hsum = H₀ + δH

#show(typeof(Hsum))