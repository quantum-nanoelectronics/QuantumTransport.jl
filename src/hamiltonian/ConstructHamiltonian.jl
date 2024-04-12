using LinearAlgebra
using SparseArrays
using Arpack
using Distributions

# For some reason, needed to add this for driver to work correctly after fixing other errors
function xyztor(p, ivec)
    ix = ivec[1]
    iy = ivec[2]
    iz = ivec[3]
    isite = ivec[4]
    δdict = Dict(0 => p["A"] * [0.0; 0.0; 0.0], #In 1 
        1 => p["A"] * [0.5; 0.5; 0.5]) #In 2
    #2 => pA*[0.0; 0.5; 0.8975-0.5], #Bi 1
    #3 => pA*[0.5; 0.0; 1.10248-0.5]) #Bi 2
    δ = δdict[isite]
    R = p["a₁"] * ix + p["a₂"] * iy + p["a₃"] * iz + δ
    return R
end

# Takes in the parameter list and the vector potential
function genH(p, A, H₀, edge_NNs, returnvals)
    ⊗(A, B) = kron(A, B)
    NormalDist = Normal(0, p["μ_disorder"])
    H_onsite = Diagonal(rand(NormalDist, p["nx"] * p["ny"] * p["nz"] * p["nsite"])) ⊗ I(p["norb"] * p["nspin"]) .+ p["μ"] * I(p["n"])
    Hᵦ = 0I(p["n"])
    H₀ = sparse(H₀ .+ H_onsite .+ Hᵦ)
    Rvals = RvalsGen(p)
    Rsurf = Vector{Float64}[]

    if (p["deviceMagnetization"] == true)
        (Bfield, Bsurf, avgB) = fieldUtils(p, A, Rsurf, Rvals, returnvals)
        Hᵦ = zeeman(map(B -> Float64.(B), Bfield), p)
    else
        Hᵦ = 0I(p["n"])
    end

    function H(k)
        if (any(isnan, k) == true)
            throw(DomainError(k, "Something broken in k vector definition! Returning NaN"))
            return
        end
        # hamiltonian describing the edges
        #build sparse as Hₑ = sparse(rows, cols, elements)
        rows = Int[]
        cols = Int[]
        elements = ComplexF64[]
        for NN in edge_NNs
            Δϕ = exp(im * k ⋅ (p["A"] * NN.N))
            for i = 1:2
                for j = 1:2
                    push!(rows, 2 * NN.b + i - 2)
                    push!(cols, 2 * NN.a + j - 2)
                    push!(elements, copy(NN.t[i, j] * Δϕ))
                end
            end
        end
        Hₑ = sparse(rows, cols, elements)
        if (size(Hₑ) == (0, 0))
            Hₑ = spzeros(ComplexF64, p["n"], p["n"])
        end
        #println("Hedges = $(size(Hₑ)), H₀ = $(size(H₀))")
        Htot = H₀ .+ Hₑ
        return dropzeros(Htot)
    end
    println("Great, SCF converged. Returning H(k).\n")
    return H
end

function zeeman(Bvals::Vector{Vector{Float64}}, p::Dict)
    # redeclared consts (need a better way)
    ħ = 1.05457E-34
    m₀ = 9.10938E-31
    ⊗(A, B) = kron(A, B)
    σ₀ = [
        1 0
        0 1
    ]
    σ₁ = [
        0 1
        1 0
    ]
    σ₂ = [
        0 -im
        im 0
    ]
    σ₃ = [
        1 0
        0 -1
    ]
    σ = Vector{Matrix}(undef, 3)
    σ[1] = σ₁;
    σ[2] = σ₂;
    σ[3] = σ₃;


    # only defined for S-like orbitals with lz = 0
    N = p["n"]
    zeeman = spzeros(ComplexF64, N, N)
    C = ħ / (2 * m₀) #sans q factor -> eV
    for ax = 1:3
        BiVals = [B[ax] for B in Bvals]
        zeeman .+= 2 * C * Diagonal(BiVals) ⊗ I(p["norb"]) ⊗ σ[ax]
    end
    return sparse(zeeman)
end

function RvalsGen(p)
    N = p["nx"] * p["ny"] * p["nz"] * p["nsite"]
    R = Vector{Vector{Float64}}(undef, N)
    for ix = 0:(p["nx"]-1)
        for iy = 0:(p["ny"]-1)
            for iz = 0:(p["nz"]-1)
                for isite = 0:(p["nsite"]-1)
                    iR = Int(1 + isite + ix * p["nsite"] + iy * p["nx"] * p["nsite"] + iz * p["ny"] * p["nx"] * p["nsite"])
                    ivec = Int.([ix, iy, iz, isite])
                    Rval = xyztor(p, ivec)
                    R[iR] = deepcopy(Rval)
                end
            end
        end
    end
    return R # needs ⊗I(pnorb)⊗I(2) for full (spinful) hilbert space
end

function fieldUtils(p, A::Function, Rsurf::Vector{Vector{Float64}}, Rvals::Vector{Vector{Float64}}, returnvals)

    ħ = 1.05457E-34
    m₀ = 9.10938E-31

    println("===== Calculating applied zeeman field properties =====")
    #checkPeriodicField(A,p) tells you if gauge field periodic with SL temporary comment
    if (p["fieldtype"] == "A")
        println("Field type = vector potential, applying onsite zeeman and peierls term")
        Bfield = Bvals(A, Rvals)
        Nsites = p["nx"] * p["ny"] * p["nz"] * p["nsite"]
        avgB = sum(Bfield) * Nsites^-1
        #avgB = [0;0;0]
        Bfield = [Bfield[i] .- avgB for i = 1:Nsites]

        Bsurf = Bvals(A, Rsurf)
        plot3Dvectors(Rvals, Bfield, [coeff * R[2] for R in Rvals], "x position (nm)", "y position (nm)", "β₂ (eV)")
        plotScatter(Rsurf, [B[1] - avgB[1] for B in Bsurf], "x position (nm)", "y position (nm)", "Bₓ (T)", "coolwarm",)
        println("ΣB field = $(sum(Bfield)) T")
        println("max B field = $(maximum(maximum.(Bfield))) T")
        return (Bfield, Bsurf, avgB)
        #return Bfield, Bsurf, avgB
        #println("Bsurf = $(Bsurf[:][:][3])")
        #Bfield = Bfield .- avgB
    elseif (p["fieldtype"] == "β")
        avgB = [0; 0; 0]
        coeff = ħ / m₀
        Bfield = A.(Rvals)
        Bsurf = A.(Rsurf)
        #=if(p["verbose"])
            println("Field type = exchange, applying onsite exchange-like term")
            println("Σβ field = $(sum(coeff*Bfield)) eV")
            println("max β field = $(coeff*maximum(maximum.(Bfield))) eV")
        end
        #plot3Dvectors(Rvals,Bfield,[coeff*B[2] for B in Bfield],"x position (nm)", "y position (nm)", "z position (nm)", "β₂ (eV)")
        if(p["plotfield"])
            colorfunc(B) = B⋅[0;1;0]
            devicefig = render2TDevice(p,Rvals,coeff*Bfield,colorfunc,10.0) 
            if("plotfield" ∈ p["returnvals"])
                push!(returnvals,devicefig)
            end
            #plotScatter(Rsurf,[coeff*B[2] for B in Bsurf], "x position (nm)", "y position (nm)", "β₂ (eV)", "coolwarm",)
        end=#
        #return (Float64.(Bfield), Float64.(Bsurf), avgB)
        return (Bfield, Bsurf, avgB)
    end
end