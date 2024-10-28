# include("Structs.jl")
function xyzElectrodeSiteToI(ElectrodeInfo::Electrode, ivec::Vector{Int})
    nx = Int(abs(ElectrodeInfo.xrange[2] - ElectrodeInfo.xrange[1]))
    ny = Int(abs(ElectrodeInfo.yrange[2] - ElectrodeInfo.yrange[1]))
    nz = Int(abs(ElectrodeInfo.zrange[2] - ElectrodeInfo.zrange[1]))
    ix = 0 * ivec[1]
    iy = ivec[2]
    iz = ivec[3]
    return ix + nx * (iy + ny * (iz)) + 1
end

function xyzElectrodeToI(p::NamedTuple, ElectrodeInfo::Electrode, ivec::Vector{Int})
    nx = Int(abs(ElectrodeInfo.xrange[2] - ElectrodeInfo.xrange[1]))
    ny = Int(abs(ElectrodeInfo.yrange[2] - ElectrodeInfo.yrange[1]))
    nz = Int(abs(ElectrodeInfo.zrange[2] - ElectrodeInfo.zrange[1]))
    ix = ivec[1]
    iy = ivec[2]
    iz = ivec[3]
    iorb = ivec[4]
    return iorb + p.norb * (ix + nx * (iy + ny * (iz))) + 1
end

function electrodeSiteToDeviceIndex(p::Dict, ElectrodeInfo::Electrode, ivecContact::Vector{Int})
    # now construct an ivec for the site in the device
    if (ElectrodeInfo.connectfrom == "-x")
        ix = 0 # the maximal site in x, edge of electrode
    else # just uhh, presuming that we will only connect in +- x. Can be changed...
        ix = (p["nx"] - 1)
    end
    iy = ivecContact[2] + ElectrodeInfo.yrange[1]
    iz = ivecContact[3] + ElectrodeInfo.zrange[1]
    return ix + p["nx"] * (iy + p["ny"] * iz) + 1
end

function changeBasis(p::Dict, ElectrodeInfo::Electrode)
    #nE = ElectrodeInfo.n*p.nsite
    #nD = p.n*p.nsite
    nE = ElectrodeInfo.n
    nD = p["nx"] * p["ny"] * p["nz"]
    #println("nE = $nE; nD = $nD")
    Psite = spzeros(nD, nE)
    nx = Int(abs(ElectrodeInfo.xrange[2] - ElectrodeInfo.xrange[1]))
    ny = Int(abs(ElectrodeInfo.yrange[2] - ElectrodeInfo.yrange[1]))
    nz = Int(abs(ElectrodeInfo.zrange[2] - ElectrodeInfo.zrange[1]))
    # only consider sites that will actually touch the slab
    ix = 0
    if (ElectrodeInfo.connectfrom == "-x")
        ix = p["nx"] - 1 # the maximal site in x, edge of electrode
    else # just uhh, presuming that we will only connect in +- x. Can be changed...
        ix = 0
    end
    for iy = 0:(ny-1)
        for iz = 0:(nz-1)
            ivec = [ix, iy, iz, 0]
            contactSiteIndex = xyzElectrodeSiteToI(ElectrodeInfo, ivec)
            deviceSiteIndex = electrodeSiteToDeviceIndex(p, ElectrodeInfo, ivec)
            #println("Device site: $deviceSiteIndex, Contact site: $contactSiteIndex")
            Psite[deviceSiteIndex, contactSiteIndex] = 1
        end
    end
    return dropzeros(Psite ⊗ I(p["nsite"] * p["norb"] * 2))
    #return sparse(Psite⊗ones(p.nsite,p.nsite)⊗I(p.norb*2))
end


function nsite(isite::Int)
    return -2 * isite + 1
end

