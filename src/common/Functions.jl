# The functions in this file are common functions that are used in multiple files.

function ⊗(A, B)
    return kron(A, B)
end

# ⊗(A, B) = kron(A, B)


# Takes in the parameters p, index vector (ix,iy,iz,isite (in unit cell), and iorb)
# returns the site-index in the full hamiltonian
function xyztoi(p, ivec, N::Vector{Int}=[0; 0; 0])
    # indexing 0 to N-1
    # in case the A₂ lattice vector != C*[0;1;0]
    diy = Int(round(p["SLa₂"][1] / p["a₁"][1])) * N[2]
    ix = mod(ivec[1], p["nx"])
    iy = mod(ivec[2] + diy, p["ny"])
    iz = mod(ivec[3], p["nz"])
    isite = ivec[4]
    iorb = ivec[5]
    return iorb + p["norb"] * isite + p["nsite"] * p["norb"] * ix + p["nsite"] * p["norb"] * p["nx"] * iy + p["nsite"] * p["norb"] * p["nx"] * p["ny"] * iz + 1
end

# Same as above, except returns the corresponding atomic position of each index vector 
# useful for calculating ∫A⋅δR peierls phase
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


×(u, v) = cross(u, v)

rot(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

function genBZ(p::Dict,nx::Int=0, ny::Int=100, nz::Int=100) # only works for cubic lattice
    # nx, ny, and nz specifically refer to # of points in IBZ
    kpoints = Vector{Float64}[]
    kindices = Vector{Int}[]
    kweights = Float64[]
    G1 = p["G"][:,1];
    G2 = p["G"][:,2];
    G3 = p["G"][:,3];
    
    if(nx != 0)
        kxs = collect(LinRange(-G1[1],G1[1],2*nx+1));
    else
        kxs = [0]
    end
    if(ny != 0)
        kys = collect(LinRange(-G2[2],G2[2],2*ny+1));
    else
        kys = [0]
    end
    if(nz != 0)
        kzs = collect(LinRange(-G3[3],G3[3],2*nz+1));
    else
        kzs = [0]
    end

    function divFixNaN(a::Int,b::Int) # for this particular instance, n/0 represents a Γ-centred sampling @ k = 0. 
        if(b==0)
            return 0
        else
            return a/b
        end
    end

    for ix = -nx:nx
        for iy = -ny:ny
            for iz = -nz:nz
                kindex = [iy + ny + 1; iz + nz + 1]
                k = divFixNaN(ix,nx)*G1/2 + divFixNaN(iy,ny)*G2/2 + divFixNaN(iz,nz)*G3
                kweight = 1
                if(abs(ix) == nx)
                    kweight *= 1/2
                end
                if(abs(iy) == ny)
                    kweight *= 1/2
                end
                if(abs(iz) == nz)
                    kweight *= 1/2
                end
                push!(kpoints,k)
                push!(kindices,kindex)
                push!(kweights,kweight)
            end
        end
    end
    ksum = sum([w for w in kweights])
    kweights = (1/ksum).*kweights
    return kpoints, kweights, kindices, kxs, kys, kzs
end
