using LinearAlgebra
using SparseArrays
include("Utilities.jl")
include("Hamiltonians.jl")

# return a vector of Σ(k) functions which return Σₖ(E) which return a sparse nsite × nsite matrix at a given energy
function genΣₖs(p::Dict, ElectrodeInfo::Vector{Electrode})
    nE = size(ElectrodeInfo)[1]
    Σks = Vector{Function}(undef, nE)
    for i = 1:nE # define a self-energy function for every electrode attached to device
        NNs = genNNs(p, ElectrodeInfo[i])
        P = changeBasis(p, ElectrodeInfo[i])
        #(Hslab, Hₗ, Hᵣ) = HcontactGen(p,NNs,ElectrodeInfo[i]) # returns in-electrode matrix, then H_inelectrode(k), H_contact(k) for transverse k
        Hs = HcontactGen(p, NNs, ElectrodeInfo[i])
        #println(Hs[1](k))
        # Thus, Hs = (Hslab(k), Hₗ(k), Hᵣ(k))
        # ∃ one contact on left, (nE-1) Contacts on right
        iCinD = Int.(sign(-i + 1.5) / 2 + 2.5) # gets the right index for the coupling hamiltonian  
        iCfromD = Int.(sign(i - 1.5) / 2 + 2.5) # gets the right index for the coupling hamiltonian  
        if (p["electrodeMaterial"] == "mtjweyl")
            # doing some jank to get the coupling matrices right
            insElectrode = deepcopy(ElectrodeInfo[i])
            insElectrode.type = "wins"
            insNNs = genNNs(p, insElectrode)
            #println("InsNNs: \n")
            #display(insNNs)
            #insHslab, insHₗ, insHᵣ = HcontactGen(p,NNs,insElectrode) # returns in-electrode matrix, then H_inelectrode(k), H_contact(k) for transverse k
            HsWMTJ = HcontactGen(p, insNNs, insElectrode)
            V = HsWMTJ[iCinD]
            βₐ = Hs[iCfromD] # away from device
            βₜ = Hs[iCinD] # towards device
        else
            V = Hs[iCinD]
            βₐ = Hs[iCfromD] # away from device
            βₜ = Hs[iCinD] # towards device
        end
        #kxes, kweights, kindices = genTetBZ(electrodeParams(p,ElectrodeInfo[1]),1000,0,0)
        Σk(k) = TΣgen(p, Hs[1](k), βₐ(k), βₜ(k), V(k), ElectrodeInfo[i], P)
        Σks[i] = Σk
    end
    return Σks
end

# implements the Sancho-Rubio method for efficient contact self-energy generation
function TΣgen(p::Dict, H::SparseMatrixCSC, βₐ::SparseMatrixCSC, βₜ::SparseMatrixCSC, V::SparseMatrixCSC, ElectrodeInfo::Electrode, P, cutoff::Float64=10^-12 * eV)
    #function Σgen(p::NamedTuple,H::SparseMatrixCSC,H_coupling::SparseMatrixCSC, Hᵥ::SparseMatrixCSC, ElectrodeInfo::Electrode, P, cutoff::Float64=10^-7*eV)
    n = ElectrodeInfo.n * p["nsite"] * p["norb"] * 2
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
        ω = (E + im * p["η"]) * I(n)
        t₀ = inv(Array(ω .- H)) * βₐ
        #t̃₀ = grInv(ω.-H)*βₜ
        t̃₀ = inv(Array(ω .- H)) * βₜ
        #t̃₀ = grInv(ω.-H)*βₜ
        #t₁ = grInv(I(n) .- t₀*t̃₀ .- t̃₀*t₀)*t₀^2
        t̃₋ = I(n)
        t̃ᵢ = Array(t̃₀)
        t₋ = I(n)
        tᵢ = Array(t₀)
        Tᵢ = zeros(ComplexF64, n, n)
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
            Π = Π * t̃₋
            #Π = deepcopy(Π*t̃₋)
            #display(Π)
            #println("")
            t̃ᵢ = inv(I(n) .- t₋ * t̃₋ .- t̃₋ * t₋) * t̃₋^2
            tᵢ = inv(I(n) .- t₋ * t̃₋ .- t̃₋ * t₋) * t₋^2
            Tᵢ .+= Π * tᵢ
            error = norm(Tᵢ .- T₋) / norm(Tᵢ)
        end
        #println("Converged, error = $error")
        effH = Array(ω .- H .- βₜ * Tᵢ)
        Gsurf = inv(effH)

        #Gsurf = grInv((E+im*p.η)*I(n) .- H .- conj(V)*Gₑ*βₐ) # guess for the green's functions in the electrodes
        #Gsurf = grInv(
        Σ_surf = V * Gsurf * V'
        #Σ_surf = V*Gsurf*V'
        #Σ_surf = (βₐ)*Gₑ*(V)'
        #Σ_surf = spzeros(ComplexF64,n,n)
        return sparse(P * Σ_surf * P')
    end
    return Σ
end

