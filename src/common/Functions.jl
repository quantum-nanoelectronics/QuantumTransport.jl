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
    X1 = p["kdict"]["X₁"];
    X2 = p["kdict"]["X₂"];
    X3 = p["kdict"]["X₃"];
    if(nx != 0)
        kxs = collect(LinRange(-X1[1],X1[1],2*nx+1));
    else
        kxs = [0]
    end
    if(ny != 0)
        kys = collect(LinRange(-X2[2],X2[2],2*ny+1));
    else
        kys = [0]
    end
    if(nz != 0)
        kzs = collect(LinRange(-X3[3],X3[3],2*nz+1));
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
                k = divFixNaN(ix,nx)*X1 + divFixNaN(iy,ny)*X2 + divFixNaN(iz,nz)*X3
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
