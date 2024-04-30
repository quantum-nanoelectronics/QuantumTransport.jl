
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


function getbands(p::Dict, klist::Vector{String}, kdict::Dict, n::Int, Hofk::Function) #takes in array of kpt strings, number of interpolation pts, H(k) function
    kpts = interpolate(n, klist, kdict)
    nk = size(kpts)[1]
    testH = Hofk([0.0; 0.0; 0.0])
    if (any(isnan, testH) == true)
        throw(DomainError(testH, "Something broken in hamiltonian definition! Returning NaN"))
        return
    end
    
    nE = size(testH)[1]
    nEig = nE
    Evals = zeros(Float64, nk, nE)
    Estates = zeros(ComplexF64, nk, nEig, nE)
    for ik in 1:nk
        k = kpts[ik]
        H = Hofk(k)
        Eofk, Estatek = eigen(Array(H))
        for iE in 1:nE
            Evals[ik, iE] = real(Eofk[iE])
        end
        for iE1 in 1:nE
            Estates[ik, :, iE1] = Estatek[:, iE1]
        end
    end

    println("Saving band structure")

    println("size of kpts: ", size(kpts))
    println("size of Evals: ", size(Evals))
    println("size of Estates: ", size(Estates))

    # TODO broken
    # save_data_formatted("bandstructure", p["path"], "bandstructure.csv", ["X", "Y"], [[0,0,0], [0,0,0], [0,0,0]]; flip_axes=true, title="Bands")
    # save_data_formatted("bandstructure", p["path"], "bandstructure.csv", ["X", "Y"], [kpts, Evals, Estates]; flip_axes=true, title="Bands")

end
