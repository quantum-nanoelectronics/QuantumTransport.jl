using LinearAlgebra
using SparseArrays
include("LowRankMatrices.jl")
using .LowRankMatrices
include("PreconditionedInverse.jl")
using .PreconditionedInverse
#import .LowRankMatrices

n = 8
H = ComplexF64.(Tridiagonal(-ones(n-1),1*ones(n),-ones(n-1)))
H = sparse(H)
H[1,n] = -1; H[n,1] = -1
#display(eigvals(H))

#H = rand(n,n)
H₀ = LRA(sparse(H), -0.2, 0.2, initialrank=4, rankstep=2)
shiftH₀ = H₀+(0.1)*I(n)
show(H₀)

δH = 0.001*rand(n,n)
display(eigvals(H+δH))
Leigvecs, Leigvals, Reigvals, Reigvecs = PreconditionedInverse.eigen((H₀, δH), H₀.Rvecs, H₀.λs)
display(Reigvals)

#=show(typeof(H₀))

δH = 1E-3*inv(H₀)
show(δH)
display(real.(Array(H₀)))
display(Array(H))
#Hsum = H₀ + δH

#show(typeof(Hsum))=#