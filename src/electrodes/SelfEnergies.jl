using LinearAlgebra
using SparseArrays

include("Structs.jl")

# return a vector of Σ(k) functions which return Σₖ(E) which return a sparse nsite × nsite matrix at a given energy
function genΣₖs(p::NamedTuple,ElectrodeInfo::Vector{Electrode})   
	nE = size(ElectrodeInfo)[1]
	Σks = Vector{Function}(undef,nE)
	for i = 1:nE # define a self-energy function for every electrode attached to device
            NNs = genNNs(p,ElectrodeInfo[i])
            P = changeBasis(p,ElectrodeInfo[i])
            #(Hslab, Hₗ, Hᵣ) = HcontactGen(p,NNs,ElectrodeInfo[i]) # returns in-electrode matrix, then H_inelectrode(k), H_contact(k) for transverse k
            Hs = HcontactGen(p,NNs,ElectrodeInfo[i]) 
            # Thus, Hs = (Hslab(k), Hₗ(k), Hᵣ(k))
            # ∃ one contact on left, (nE-1) Contacts on right
            iCinD = Int.(sign(-i+1.5)/2+2.5) # gets the right index for the coupling hamiltonian  
            iCfromD = Int.(sign(i-1.5)/2+2.5) # gets the right index for the coupling hamiltonian  
            if(p.electrodeMaterial=="mtjweyl")
                # doing some jank to get the coupling matrices right
                insElectrode = deepcopy(ElectrodeInfo[i])
                insElectrode.type = "wins"
                insNNs = genNNs(p,insElectrode)
                #println("InsNNs: \n")
                #display(insNNs)
                #insHslab, insHₗ, insHᵣ = HcontactGen(p,NNs,insElectrode) # returns in-electrode matrix, then H_inelectrode(k), H_contact(k) for transverse k
                HsWMTJ = HcontactGen(p,insNNs,insElectrode)
                V  = HsWMTJ[iCinD]
                βₐ = Hs[iCfromD] # away from device
                βₜ = Hs[iCinD] # towards device
            else
                V  = Hs[iCinD]
                βₐ = Hs[iCfromD] # away from device
                βₜ = Hs[iCinD] # towards device
            end
            #Σk(k) = Σgen(p,Hs[1](k),Hᵥₑ(k),Hᵥₘ(k),ElectrodeInfo[i],P)
            kxes, kweights, kindices = genTetBZ(electrodeParams(p,ElectrodeInfo[1]),1000,0,0)
            Σk(k) = TΣgen(p,Hs[1](k),βₐ(k),βₜ(k),V(k),ElectrodeInfo[i],P)
            #Σk(k) = Σgen(p,Hs[1](k),βₐ(k),βₜ(k),V(k),ElectrodeInfo[i],P)
            #Σk(k) = ∫Σgen(p,Hs[4], Hs[1](k), βₐ(k), βₜ(k), V(k),ElectrodeInfo[i],P,k,kxes,kweights)
            #Σk(k) = Σgen(p,Hs[1](k),Hs[2](k).+Hs[3](k),Hᵥₘ(k),ElectrodeInfo[i],P)
            Σks[i] = Σk
	end
	return Σks
end

#=function genΣₖs(p::NamedTuple,ElectrodeInfo::Vector{Electrode})   
	nE = size(ElectrodeInfo)[1]
	Σks = Vector{Function}(undef,nE)
	for i = 1:nE # define a self-energy function for every electrode attached to device
            if(p.electrodeMaterial=="mtjweyl")
                # for the normal part of the weyl electrodes
                NNs = genNNs(p,ElectrodeInfo[i])
		P = changeBasis(p,ElectrodeInfo[i])
		Hslab, Hₗ, Hᵣ = HcontactGen(p,NNs,ElectrodeInfo[i]) # returns in-electrode matrix, then H_inelectrode(k), H_contact(k) for transverse k
                # doing some jank to get the coupling matrices right
                insElectrode = deepcopy(ElectrodeInfo[i])
                insElectrode.type = "wins"
                insNNs = genNNs(p,insElectrode)
		insHslab, insHₗ, insHᵣ = HcontactGen(p,NNs,insElectrode) # returns in-electrode matrix, then H_inelectrode(k), H_contact(k) for transverse k
                connectDict = Dict{String,Function}("-x"=>insHᵣ,"+x"=>insHₗ)
                Hᵥ = connectDict[ElectrodeInfo[i].connectfrom]
                Σk(k) = Σgen(p,Hslab(k),Hₗ(k).+Hᵣ(k),Hᵥ(k),ElectrodeInfo[i],P)
		Σks[i] = Σk
            else
                NNs = genNNs(p,ElectrodeInfo[i])
		P = changeBasis(p,ElectrodeInfo[i])
		Hslab, Hₗ, Hᵣ = HcontactGen(p,NNs,ElectrodeInfo[i]) # returns in-electrode matrix, then H_inelectrode(k), H_contact(k) for transverse k
                #connectDict = Dict{String,Function}("-x"=>Hᵣ,"+x"=>Hₗ)
                if(ElectrodeInfo[i].connectfrom=="-x")
                    Hᵥ = Hᵣ
                    Σk = k -> Σgen(p,Hslab(k),Hₗ(k).+Hᵣ(k),Hᵥ(k),ElectrodeInfo[i],P)
                    #Σk(k) = Σgen(p,Hslab(k),Hₗ(k).+Hᵣ(k),Hᵥ(k),ElectrodeInfo[i],P)
                    Σks[i] = deepcopy(Σk)
                    #Σk(k) = Σgen(p,Hslab(k),Hₗ(k).+Hᵣ(k),Hᵥ(k),ElectrodeInfo[i],P)
		    #Σks[i] = Σk
                else
                    Hᵥ = Hₗ
                    Σk = k -> Σgen(p,Hslabk),Hₗ(k).+Hᵣ(k),Hᵥ(k),ElectrodeInfo[i],P)
                    #Σk(k) = Σgen(p,Hslab(k),Hₗ(k).+Hᵣ(k),Hᵥ(k),ElectrodeInfo[i],P)
                    Σks[i] = deepcopy(Σk)
                end
                #Hᵥ = connectDict[ElectrodeInfo[i].connectfrom]
            end
            #Σk(k) = Σgen(p,Hslab(k),Hₗ(k).+Hᵣ(k),Hᵥ(k),ElectrodeInfo[i],P)
            #Σk(k) = Σgen(p,Hslab(k),Hₗ(k),Hᵣ(k),ElectrodeInfo[i],P)
	end
	return Σks
end=#

function TΣgen(p::NamedTuple,H::SparseMatrixCSC,βₐ::SparseMatrixCSC, βₜ::SparseMatrixCSC, V::SparseMatrixCSC, ElectrodeInfo::Electrode, P, cutoff::Float64=10^-12*eV)
#function Σgen(p::NamedTuple,H::SparseMatrixCSC,H_coupling::SparseMatrixCSC, Hᵥ::SparseMatrixCSC, ElectrodeInfo::Electrode, P, cutoff::Float64=10^-7*eV)
    n = ElectrodeInfo.n*p.nsite*p.norb*2
    # so H needs to be instantiated and called outside of the loop
    #H_coupling = Hₗ .+ Hᵣ# couples a layer to the infinite on both sides
    #BLAS.set_num_threads(1) # disable linalg multithreading and parallelize over k instead
    function Σ(E::Float64)
        #Gₑ = grInv((E+im*p.η)*I(n) .- H) # guess for the green's functions in the electrodes
        #Σ_guess = H_coupling*grInv((E+im*p.η)*I(n) .- H .- 0.1*I(n))*Hᵥₘ'
        #Σ_guess = H_coupling*grInv((E+im*p.η)*I(n) .- H .- 0.1*I(n))*H_coupling'
        # converge the self energy 
        error = 1
        # using transfer matrix method described in 
        ω = (E + im*p.η)*I(n)
        t₀ = inv(Array(ω.-H))*βₐ
        #t̃₀ = grInv(ω.-H)*βₜ
        t̃₀ = inv(Array(ω.-H))*βₜ
        #t̃₀ = grInv(ω.-H)*βₜ
        #t₁ = grInv(I(n) .- t₀*t̃₀ .- t̃₀*t₀)*t₀^2
        t̃₋ = I(n); t̃ᵢ = Array(t̃₀)
        t₋ = I(n); tᵢ = Array(t₀)
        Tᵢ = zeros(ComplexF64,n,n)
        Tᵢ .+= t₀
		Π = I(n)
		while error > cutoff
			#println("error = $error, size(T) = $(size(Tᵢ))")
			T₋ = deepcopy(Tᵢ)
			t̃₋ = t̃ᵢ
			t₋ = tᵢ
			#T₋ = deepcopy(Tᵢ)
			#t̃₋ = deepcopy(t̃ᵢ)
			#t₋ = deepcopy(tᵢ)
			Π = Π*t̃₋
			#Π = deepcopy(Π*t̃₋)
			#display(Π)
			#println("")
        	t̃ᵢ = inv(I(n) .- t₋*t̃₋ .- t̃₋*t₋)*t̃₋^2
        	tᵢ = inv(I(n) .- t₋*t̃₋ .- t̃₋*t₋)*t₋^2
                Tᵢ .+= Π*tᵢ
			error = norm(Tᵢ .- T₋)/norm(Tᵢ)	
		end
                #println("Converged, error = $error")
                effH = Array(ω .- H .- βₜ*Tᵢ)
		Gsurf = inv(effH)

		#Gsurf = grInv((E+im*p.η)*I(n) .- H .- conj(V)*Gₑ*βₐ) # guess for the green's functions in the electrodes
        #Gsurf = grInv(
        Σ_surf = V*Gsurf*V'
        #Σ_surf = V*Gsurf*V'
        #Σ_surf = (βₐ)*Gₑ*(V)'
        #Σ_surf = spzeros(ComplexF64,n,n)
        return sparse(P*Σ_surf*P')
    end
    return Σ
end
#function Σgen(p::NamedTuple,H::Matrix,Hₗ::Matrix, Hᵣ::Matrix, ElectrodeInfo::Electrode, cutoff::Float64=10^-7*eV)
function Σgen(p::NamedTuple,H::SparseMatrixCSC,βₐ::SparseMatrixCSC, βₜ::SparseMatrixCSC, V::SparseMatrixCSC, ElectrodeInfo::Electrode, P, cutoff::Float64=10^-13*eV)
#function Σgen(p::NamedTuple,H::SparseMatrixCSC,H_coupling::SparseMatrixCSC, Hᵥ::SparseMatrixCSC, ElectrodeInfo::Electrode, P, cutoff::Float64=10^-7*eV)
    n = ElectrodeInfo.n*p.nsite*p.norb*2
    # so H needs to be instantiated and called outside of the loop
    #H_coupling = Hₗ .+ Hᵣ# couples a layer to the infinite on both sides
    #BLAS.set_num_threads(1) # disable linalg multithreading and parallelize over k instead
    function Σ(E::Float64)
        #Gₑ = grInv((E+im*p.η)*I(n) .- H .- Hᵥₐ'*I(n)*Hᵥₘ) # guess for the green's functions in the electrodes
        Gₑ = grInv((E+im*p.η)*I(n) .- H) # guess for the green's functions in the electrodes
        #Σ_guess = H_coupling*grInv((E+im*p.η)*I(n) .- H .- 0.1*I(n))*Hᵥₘ'
        #Σ_guess = H_coupling*grInv((E+im*p.η)*I(n) .- H .- 0.1*I(n))*H_coupling'
        # converge the self energy 
        error = 1
        Σ = βₜ*grInv((E+im*p.η)*I(n) .- H)*βₐ # guess for the green's functions in the electrodes
        #β₁ = βₜ; β₂ = βₐ
        #β = βₜ .+ βₐ
        while error > cutoff
            #println("Gₑʳ convergence error loop: $error")
            Σ₀ = copy(Σ)
            Σ = βₜ*grInv((E+im*p.η)*I(n) .- H .- Σ₀)*βₐ # guess for the green's functions in the electrodes
            #Σ = βₜ*Gₑ0*βₐ
            #Gₑ = grInv((E+im*p.η)*I(n) .- H .- Σ_contact) # guess for the green's functions in the electrodes
            #Gₑ = grInv((E+im*p.η)*I(n) .- H .- conj(βₜ)*Gₑ0*βₐ) # guess for the green's functions in the electrodes
            #Gₑ = grInv((E+im*p.η)*I(n) .- H .- βₜ'*Gₑ0*βₐ) # guess for the green's functions in the electrodes
            #Gₑ = grInv((E+im*p.η)*I(n) .- H .- β'*Gₑ0*β) # guess for the green's functions in the electrodes
            #Gₑ = grInv((E+im*p.η)*I(n) .- H .- β'*Gₑ0*β) # guess for the green's functions in the electrodes
            #Gₑ = grInv((E+im*p.η)*I(n) .- H .- β'*Gₑ0*β) # guess for the green's functions in the electrodes
            error =  norm(Σ.-Σ₀)/norm(Σ₀)
        end
        #println("\nΣ = ")
        # loop to SCF for surface?
        #Σ_surf = (V.+βₐ)*Gₑ*(V.+βₐ)'
        #Gsurf = grInv((E+im*p.η)*I(n) .- H .- V*Gₑ*βₐ) # guess for the green's functions in the electrodes
        #Gsurf = grInv((E+im*p.η)*I(n) .- H .- conj(V)*Gₑ*βₐ) # guess for the green's functions in the electrodes
        #Gsurf = grInv(
        Σ_surf = V*Gₑ*βₐ
        #Σ_surf = (βₐ)*Gₑ*(V)'
        #Σ_surf = βₐ*Gₑ*βₜ'
        #Σ_surf = βₜ*Gₑ*βₐ'
        #Σ_surf = V*Gₑ*βₐ'
        #Σ_surf = Hᵥ*grInv((E+im*p.η)*I(n) .- H .- Σ_guess)*Hᵥ'
        #Σ_surf = Hᵥ'*grInv((E+im*p.η)*I(n) .- H .- Σ_guess)*Hᵥ
        #Σ_surf = spzeros(ComplexF64,n,n)
        return P*Σ_surf*P'
        #for i = 1:40
        #=while error > cutoff*n
            #println("Σ convergence error loop: $error")
            Σ_guess0 = deepcopy(Σ_guess)
            Σ_guess = H_coupling*grInv((E+im*p.η)*I(n) .- H .- Σ_guess0)*H_coupling'
            #Σ_guess = H_coupling'*grInv((E+im*p.η)*I(n) .- H .- Σ_guess0)*H_coupling
            error =  norm(Σ_guess.-Σ_guess0)
        end
        #println("\nΣ = ")
        # loop to SCF for surface?
        Σ_surf = Hᵥ*grInv((E+im*p.η)*I(n) .- H .- Σ_guess)*Hᵥ'
        #Σ_surf = Hᵥ'*grInv((E+im*p.η)*I(n) .- H .- Σ_guess)*Hᵥ
        return P*Σ_surf*P'=#
        #=while error > cutoff*n
            #println("Σ convergence error loop: $error")
            Σ_surf0 = deepcopy(Σ_surf)
            Σ_surf = Hᵥ*grInv((E+im*p.η)*I(n) .- H .- Σ_guess .-Σ_surf0)*Hᵥ'
            error =  norm(Σ_surf.-Σ_surf0)
        end
        #display(Σ_guess)
        #Σ_surf = Hᵥ*grInv((E+im*p.η)*I(n) .- H .- Σ_guess)*Hᵥ'
        return P*Σ_surf*P'=#
        #=if(ElectrodeInfo.connectfrom=="-x")
                Σ_surf = Hᵣ*grInv((E+im*p.η)*I(n) .- H .- Σ_guess)*Hᵣ'
                return P*Σ_surf*P'
        else
                if(ElectrodeInfo.connectfrom != "+x")
                        println("Something is very broken, check that your electrodes are in ±x")
                end
                Σ_surf = Hₗ*grInv((E+im*p.η)*I(n) .- H .- Σ_guess)*Hₗ'
                #show(size(Σ_surf))
                #show(size(P))
                return P*Σ_surf*P'
        end=#
    end
    return Σ
end


function ∫Σgen(p::NamedTuple, H::Function, Hslab::SparseMatrixCSC, βₐ::SparseMatrixCSC, βₜ::SparseMatrixCSC, V::SparseMatrixCSC, ElectrodeInfo::Electrode, P::SparseMatrixCSC, k::Vector{Float64}, kxes::Vector{Vector{Float64}}, kweights::Vector{Float64}, cutoff::Float64=10^-5*eV)
    n = ElectrodeInfo.n*p.nsite*p.norb*2
    kvals = [[kₓ[1],k[2],k[3]] for kₓ in kxes]
    function Σ(E::Float64)
        Gₑ = spzeros(ComplexF64,n,n)
        Gₑs = map(k->grInv((E+im*p.η)*I(n) .- H(k)),kvals)
        for i in eachindex(kweights)
                Gₑ .+= kweights[i]*Gₑs[i]
        end
        Gsurf = grInv((E+im*p.η)*I(n) .- Hslab .- βₜ*Gₑ*βₐ)
        Σ_surf = V*Gₑ*V'
        return P*Σ_surf*P'
    end
    return Σ
end

