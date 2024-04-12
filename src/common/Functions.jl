# The functions in this file are common functions that are used in multiple files.

f(x) = x^2

function ⊗(A, B)
    return kron(A, B)
end

# ⊗(A, B) = kron(A, B)


# ×(u, v) = cross(u, v)

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
