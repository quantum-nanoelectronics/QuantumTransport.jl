using LinearAlgebra
using SparseArrays
using Arpack
using Distributions

# Takes in the parameter list and the vector potential
function genH(p, H₀, edge_NNs)
     NormalDist = Normal(0,p.μ_disorder)
     H_onsite = Diagonal(rand(NormalDist,p.n*p.nsite))⊗I(p.norb*2) .+ p.μ*I(p.n*p.nsite*p.norb*2)
     Hᵦ = 0I(p.n*p.nsite*p.norb*2)
	ntot = p.n*p.nsite*p.norb*2
	H₀ = sparse(H₀ .+ H_onsite .+ Hᵦ)
	function H(k)
		if(any(isnan,k)==true)
               throw(DomainError(k, "Something broken in k vector definition! Returning NaN"))
               return
		end
		# hamiltonian describing the edges
          #build sparse as Hₑ = sparse(rows, cols, elements)
          rows = Int[]; cols = Int[]; elements = ComplexF64[];
          for NN in edge_NNs
               Δϕ = exp(im*k⋅(p.A*NN.N))
               for i = 1:2
                    for j = 1:2
                         push!(rows,2*NN.b+i-2);
                         push!(cols,2*NN.a+j-2);
                         push!(elements,copy(NN.t[i,j]*Δϕ))
                    end
               end
		end
		Hₑ = sparse(rows,cols,elements)
		if(size(Hₑ) == (0,0)) Hₑ = spzeros(ComplexF64,ntot,ntot) end 
		#println("Hedges = $(size(Hₑ)), H₀ = $(size(H₀))")
		Htot = H₀ .+ Hₑ
		return dropzeros(Htot)
	end
	println("Great, SCF converged. Returning H(k).\n")
	return H 
end