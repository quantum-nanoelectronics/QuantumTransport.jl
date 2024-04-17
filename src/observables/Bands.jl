
function interpolate(n::Int, klist::Vector{String}, kdict::Dict)
    highSym = (1) * map(k -> kdict[k], klist)
    kpts = Any[]
    d = 1 / n
    for k1 in 1:(size(highSym)[1]-1)
        k2 = k1 + 1
        for f = 0:d:(1-d) #fraction of the 2nd kpt
            k = (1 - f) * highSym[k1] + f * highSym[k2]
            push!(kpts, k)
        end
    end
    push!(kpts, highSym[end])
    return kpts
end

function getbands(klist::Vector{String}, kdict::Dict, n::Int, a::Float64, Hofk::Function, arpack=false) #takes in array of kpt strings, number of interpolation pts, H(k) function
    kpts = interpolate(n, klist, kdict)
    nk = size(kpts)[1]
    testH = Hofk([0.0; 0.0; 0.0])
    if (any(isnan, testH) == true)
        throw(DomainError(testH, "Something broken in hamiltonian definition! Returning NaN"))
        return
    end
    #initialize the ban array
    #λ_test, evecs_test = eig(testH, 1E-12)
    #arpack = true # small hamiltonian, few bands
    #=if (arpack != false)
        maxiter = 8000
        #nE = 128
        nE = arpack
        #nE = 6*Int(floor(log2(size(testH)[1])))
        nEig = size(testH)[1]
        if (nE < size(testH)[1])
            println("Heads up! $nE / $nEig eigvls are being calculated")
        end
        λ_test, evecs_test = eigs(testH, nev=nE, which=:SM, maxiter=maxiter)
        #nE = size(λ_test)[1]
        ndim = size(evecs_test)[2]
    else=#
    nE = size(testH)[1]
    nEig = nE
    #end
    Evals = zeros(Float64, nk, nE)
    Estates = zeros(ComplexF64, nk, nEig, nE)
    for ik in 1:nk
        k = kpts[ik]
        H = Hofk(k)
        #Eofk, Estatek = eigs(Hermitian(H))
        #Eofk, Estatek = eigen(H)
        #print("$(round.(k,sigdigits=3)).. ")
        #print("H = $H\n")
        #if(norm(H) < 0.01 || k⋅k≈0)
        #	Estatek = (1/√(nEig))*ones(nEig,nE); Eofk = zeros(nE)
        #=if (arpack != false && !(k ⋅ k ≈ -1))
            #println("H(k) diagonalization time:")
            Eofk, Estatek = eigs(H, nev=nE, which=:SM, maxiter=maxiter)
            #@time Eofk, Estatek = eigs(H,nev=nE, which=:SM, maxiter=maxiter)
            #println("Done timing.")
        else
            Eofk, Estatek = eigen(Array(H))
        end =#
        #show(size(Estatek))
        Eofk, Estatek = eigen(Array(H))
        #Eofk = eigvals(H)
        for iE in 1:nE
            Evals[ik, iE] = real(Eofk[iE])
        end
        #Estatesk = eigvecs(H)
        for iE1 in 1:nE
            #for iE1 in 1:nEig; for iE2 in 1:nEig
            Estates[ik, :, iE1] = Estatek[:, iE1] # pre-arpack
            #Estates[ik,iE1,iE2] = Estatek[iE2,iE2]
        end
    end
    return kpts, Evals, Estates
end
