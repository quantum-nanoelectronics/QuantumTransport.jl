include("Structs.jl")

function xyztoi(p,ivec, N::Vector{Int} = [0;0;0]) 
	# indexing 0 to N-1
	# in case the A₂ lattice vector != C*[0;1;0]
	diy = Int(round(p.SLa₂[1]/p.a₁[1]))*N[2]
	#ix = ivec[1]; iy = ivec[2]; iz = ivec[3]; isite = ivec[4]; iorb = ivec[5]
	ix = mod(ivec[1],p.nx); iy = mod(ivec[2] + diy,p.ny); iz = mod(ivec[3],p.nz); isite = ivec[4]; iorb = ivec[5]
	#ix = mod(ivec[1],p.nx); iy = mod(ivec[2],p.ny); iz = mod(ivec[3],p.nz); isite = ivec[4]; iorb = ivec[5]
	return iorb + p.norb*isite + p.nsite*p.norb*ix + p.nsite*p.norb*p.nx*iy + p.nsite*p.norb*p.nx*p.ny*iz + 1
end

# Same as above, except returns the corresponding atomic position of each index vector 
# useful for calculating ∫A⋅δR peierls phase
function xyztor(p,ivec)
	ix = ivec[1]; iy = ivec[2]; iz = ivec[3]; isite = ivec[4];
	δdict = Dict(0 => p.A*[0.0; 0.0; 0.0], #In 1 
		     1 => p.A*[0.5; 0.5; 0.5]) #In 2
		     #2 => p.A*[0.0; 0.5; 0.8975-0.5], #Bi 1
		     #3 => p.A*[0.5; 0.0; 1.10248-0.5]) #Bi 2
	δ = δdict[isite]
	R = p.a₁*ix + p.a₂*iy + p.a₃*iz + δ
        #println("ivec = $ivec, Rpos = $R")
        return R
end

function RvalsGen(p::NamedTuple,ei::Electrode)
	N = ei.n*p.nsite
	R = Vector{Vector{Float64}}(undef,N)
	ix = -1
	if(ei.connectfrom=="-x")
		ix = -1
	elseif(ei.connectfrom=="+x")
		ix = p.nx+1
	else
		return 0 #something is busted
	end
	# offset to fix the 
	iRoffset = Int(0 + 0 + ix*p.nsite + ei.yrange[1]*p.nx*p.nsite + ei.zrange[1]*p.ny*p.nx*p.nsite)
        nx = 1
        for iy = ei.yrange[1]:(ei.yrange[2]-1)
            for iz = ei.zrange[1]:(ei.zrange[2]-1)
                for isite = 0:(p.nsite-1)
                    iR = Int(1 + isite + p.nsite*(ix-ix + nx*((iy-ei.yrange[1]) + (iz-ei.zrange[1])*p.ny)))
                    #iR = Int(1 - iRoffset + isite + ix*p.nsite + iy*p.nx*p.nsite + iz*p.ny*p.nx*p.nsite)
                    #println("$ix $iy $iz $isite iR $iR")
                    ivec = Int.([ix,iy,iz,isite])
                    Rval = xyztor(p,ivec)
                    R[iR] = deepcopy(Rval)
                end
            end
	end
	return R # needs ⊗I(p.norb)⊗I(2) for full (spinful) hilbert space
end




function xyzElectrodeSiteToI(ElectrodeInfo::Electrode,ivec::Vector{Int})
	nx = Int(abs(ElectrodeInfo.xrange[2] - ElectrodeInfo.xrange[1]))
	ny = Int(abs(ElectrodeInfo.yrange[2] - ElectrodeInfo.yrange[1]))
	nz = Int(abs(ElectrodeInfo.zrange[2] - ElectrodeInfo.zrange[1]))
	ix = 0*ivec[1]; iy = ivec[2]; iz = ivec[3]
	return  ix + nx*(iy + ny*(iz)) + 1
end

function xyzElectrodeToI(p::NamedTuple, ElectrodeInfo::Electrode,ivec::Vector{Int})
	nx = Int(abs(ElectrodeInfo.xrange[2] - ElectrodeInfo.xrange[1]))
	ny = Int(abs(ElectrodeInfo.yrange[2] - ElectrodeInfo.yrange[1]))
	nz = Int(abs(ElectrodeInfo.zrange[2] - ElectrodeInfo.zrange[1]))
	ix = ivec[1]; iy = ivec[2]; iz = ivec[3]; iorb = ivec[4]
	return  iorb + p.norb*(ix + nx*(iy + ny*(iz))) + 1
end

function electrodeSiteToDeviceIndex(p::NamedTuple, ElectrodeInfo::Electrode,ivecContact::Vector{Int})
	# now construct an ivec for the site in the device
	if(ElectrodeInfo.connectfrom=="-x")
		ix = 0 # the maximal site in x, edge of electrode
	else # just uhh, presuming that we will only connect in +- x. Can be changed...
		ix = (p.nx-1)
	end
	iy = ivecContact[2]+ElectrodeInfo.yrange[1]
	iz = ivecContact[3]+ElectrodeInfo.zrange[1]
	return ix + p.nx*(iy+p.ny*iz) + 1
end

function changeBasis(p::NamedTuple,ElectrodeInfo::Electrode)
	#nE = ElectrodeInfo.n*p.nsite
	#nD = p.n*p.nsite
	nE = ElectrodeInfo.n
	nD = p.n
        #println("nE = $nE; nD = $nD")
	Psite = spzeros(nD,nE)
	nx = Int(abs(ElectrodeInfo.xrange[2] - ElectrodeInfo.xrange[1]))
	ny = Int(abs(ElectrodeInfo.yrange[2] - ElectrodeInfo.yrange[1]))
	nz = Int(abs(ElectrodeInfo.zrange[2] - ElectrodeInfo.zrange[1]))
	# only consider sites that will actually touch the slab
	ix = 0
	if(ElectrodeInfo.connectfrom=="-x")
		ix = p.nx-1 # the maximal site in x, edge of electrode
	else # just uhh, presuming that we will only connect in +- x. Can be changed...
		ix = 0
	end
	for iy = 0:(ny-1)
		for iz = 0:(nz-1)
			ivec = [ix,iy,iz,0]
			contactSiteIndex = xyzElectrodeSiteToI(ElectrodeInfo,ivec)
			deviceSiteIndex = electrodeSiteToDeviceIndex(p,ElectrodeInfo,ivec)
			#println("Device site: $deviceSiteIndex, Contact site: $contactSiteIndex")
			Psite[deviceSiteIndex,contactSiteIndex] = 1
		end
	end
        return dropzeros(Psite⊗I(p.nsite*p.norb*2))
        #return sparse(Psite⊗ones(p.nsite,p.nsite)⊗I(p.norb*2))
end
			



function electrodeParams(p::NamedTuple,ElectrodeInfo::Electrode)
	nx = Int(abs(ElectrodeInfo.xrange[2] - ElectrodeInfo.xrange[1]))
	ny = Int(abs(ElectrodeInfo.yrange[2] - ElectrodeInfo.yrange[1]))
	nz = Int(abs(ElectrodeInfo.zrange[2] - ElectrodeInfo.zrange[1]))
	ep = (nx = nx, ny = ny, nz = nz, n = nz*nx*ny)
	return merge(p,ep)
end

rot(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

function nsite(isite::Int)
	return -2*isite + 1
end



