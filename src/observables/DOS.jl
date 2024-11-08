#function DOS(genA::Function,kvals::Vector{Vector{Float64}},kpts::Vector{Float64},E_samples::Vector{Float64})
#    DOS_samples = pmap(E->tr(genA(E))/π,E_samples)
#
#    return DOS_samples
#end


function genDOS(p::Dict, F::Function, type::Symbol, η::Float64)
    if type == :Gʳ
        Gᴿ = F    
        function DOS_Gᴿ(E::Float64)
            return (-1/π) * trace(imag.(Gᴿ))
        end
        return DOS_Gᴿ
    elseif type == :A
        A = F
        function DOS_A(E::Float64)
            return trace(A) / (2*π)
        end
        return DOS_A
    elseif type == :H
        H = F # we've got the hamiltonian as a function of K
        kpoints, dk = genBZ(p["G"],p["nkvec"])
        energies = Real[]
        for k ∈ kpoints
            append!(energies, real.(eigvals(H(k))))
        end
        # implement later
        function DOS_H(E::Float64)
            # implementation for DOS using H
            return -(dk/π)*imag(sum((-energies.+(E-im*η)).^(-1)))
        end
        return DOS_H
    end
end

# This function needed to be re-written above.
# WARNING: Method definition DOS(Float64) in module Observables at e:\WindowsSSD\ECE 464K Senior Design\QuantumTransport\src\observables\DOS.jl:10 overwritten at e:\WindowsSSD\ECE 464K Senior Design\QuantumTransport\src\observables\DOS.jl:16.
#   ** incremental compilation may be fatally broken for this module **
# function genDOS(p::Dict, F::Function, type::Symbol, η::Float64)
#     if type == :Gʳ
#         Gᴿ = F    
#         function DOS(E::Float64)
#             return (-1/π)*trace(imag.(Gᴿ))
#         end
#         return DOS
#     elseif type == :A
#         A = F
#         function DOS(E::Float64)
#             return trace(A)/(2*π)
#         end
#         return DOS
#     elseif type == :H
#         H = F # we've got the hamiltonian as a function of K
#         # implement later
#     end
# end

function genDOS(p::Dict, H::SparseMatrixCSC, η::Float64, neigs::Int, centerE::Float64)
    eigvals = eigs(H, nev = neigs, sigma=centerE)
    function DOS_fullH(E::Float64)
        return (-1/π)*sum(imag.((eigvals.+(η-E)).^-1))
    end
    return DOS_fullH
end

function genDOS(p::Dict, H::Matrix{ComplexF64}, η::Float64)
    eigvals = eigvals(H)
    function DOS_denseH(E::Float64)
        return (-1/π)*sum(imag.((eigvals.+(η-E)).^-1))
    end
    return DOS_denseH
end



function genBZ(G::Matrix, nkvec::Vector{Int})
    kpoints = Vector{Float64}[]
    nG1 = nkvec[1]; nG2 = nkvec2[2]; nG3 = nkvec3[3];
    dk = det(G)/(nG1*nG2*nG3)
    for ix ∈ 1:nG1
        for iy ∈ 1:nG2
            for iz in 1:nG3
                k = G*[ix/nG1; iy/nG2; iz/nG3]
                push!(kpoints,k)
            end
        end
    end
    return kpoints, weight
end

#=function DOS(genA::Function,kgrid::Vector{Vector{Float64}},kweights::Vector{Float64},Evals::Vector{Float64},parallelk::Bool=true)
	nE = size(Evals)
	#nkz = maximum([kindex[2] for kindex in kindices])
	#nky = maximum([kindex[1] for kindex in kindices])
	#special slice to map BZ
	nk = size(kweights)[1]
	#Eslice = findnearest(Evals,Eslice) #make it so that we do not have to do a whole nother k loop
	#TmapList = zeros(nk)
	DOS = zeros(nE)
	if(parallelk)
		#BLAS.set_num_threads(1) # disable linalg multithreading and parallelize over k instead
		iter = ProgressBar(1:nk)
		knum = shuffle([i for  i = 1:nk])
            #for ik in iter
            Threads.@threads for ik in iter
                i = knum[ik] # this is to shuffle the kpt allocations so no processor gets a dense section of grid
                k = kgrid[i]
                w = kweights[i]
                Aₖ = genA(k)
                for iE in eachindex(Evals)
                    E = Evals[iE]
                    Dₖ = deepcopy(tr(Aₖ(E)))/π
                    #TofE[iE] += real(Dₖ*w)
                    DOS[iE] += real(im*Dₖ*w)
                end
		end
	else
		#BLAS.set_num_threads(1) # disable linalg multithreading and parallelize over k instead
		knum = shuffle([i for  i = 1:nk])
		for ik = 1:nk
            i = knum[ik] # this is to shuffle the kpt allocations so no processor gets a dense section of grid
            k = kgrid[i]
            #kindex = kindices[i]
            w = kweights[i]
            Aₖ = genA(k)
            iter = ProgressBar(1:size(Evals)[1])
            for iE in iter
            #Threads.@threads for iE in iter
                E = Evals[iE]
                Dₖ = deepcopy(tr(Aₖ(E)))/π
                #DOS[iE] += real(Dₖ*w)
                DOS[iE] += real(im*Dₖ*w)
                #=if(E≈Eslice)
                        Tmap[kindex[1],kindex[2]] = real(T)
                        TmapList[i] = real(T)
                end=#
            end
		end
	end
	return DOS
end

function sitePDOS(p::NamedTuple, genGʳ::Function, Os, E::Float64=0.1*eV)
    function DOS(k::Vector{Float64})
        totToSite = sparse(I(p.nx*p.ny*p.nz)⊗(ones(p.norb*2)')) 
        Gʳ = genGʳ(k)(E);
        #SiteGᴿ = totToSite*(diag(Gᴿ))
        #display(Gᴿ); println(""); display(SiteGᴿ); println(""); display(totToSite)
        
        #DOS = (-1/π)*imag.(SiteGᴿ)
        return [Array((-1/π)*imag.(totToSite*diag(O*Gʳ))) for O in Os]
    end
    return DOS
end

function siteDOS(p::NamedTuple, genGᴿ::Function, E::Float64=0.1*eV)
    function DOS(k::Vector{Float64})
        totToSite = sparse(I(p.nx*p.ny*p.nz)⊗(ones(p.norb*2)')) 
		Gᴿ = genGᴿ(k)(E);
		SiteGᴿ = totToSite*(diag(Gᴿ))
		#display(Gᴿ); println(""); display(SiteGᴿ); println(""); display(totToSite)
		DOS = (-1/π)*imag.(SiteGᴿ)
            return Array(DOS)
    end
    return DOS
end


			
function DOS(E::Float64,A::Function,Q::Vector = I(size(A(0))[1]))
	return (1/(2*π))*tr(Q*A)
end
=#