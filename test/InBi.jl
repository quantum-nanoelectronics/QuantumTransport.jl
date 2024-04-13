module InBi
export params, kdictGen, genSL



# constants

ħ = 1 #hbar
h = ħ * 2*π
m₀ = 1 #electron mass
q = 1 #electron charge
ϵ₀ = 1 #permittivity
r₀ = 1 #bohr radii
Ry = 1/2 #rydberg
au = 1 #hartree energy
Å = 1.8897 * r₀	#angstrom, in hartree units
nm = 10 * Å 		#angstrom, in hartree units
qHo = 67;
eV = au/27.2
ε = 0.0*eV; ε₂ = 0.0*eV #onsite energy for In, Bi

# electronic properties
t = 1*eV
t₁ = 1*eV; # In-Bi px-px,py-py
t₂ = 1*eV; # In-In 1st NN hopping
t₃ = 1*eV; #Bi-Bi 1st NN hopping
t₄ = 0 # obsolete
t₅ = 0.17*eV # In-In 2nd NN hopping
t₆ = 0.04*eV;# Bi-Bi 2nd NN hopping
t₇ = 0.0*eV;# In-Bi further hoppings
t₈ = -0.1*eV;# In-In 2nd nn vertical hoppings
t₉ = 0.0*eV;# In-In px -px, py-py hoppings
ε₁ = 1.9*eV; ε₂ = -1.2*eV #onsite energy for In, Bi

# structural properties
#a = 4.8*Å; c = 4.8*Å
a = 1; c = 1
A = [a 0 0; 0 a 0; 0 0 c]
a₁ = A[:,1]; a₂ = A[:,2]; a₃ = A[:,3]
B = transpose(2*π*inv(A))
b₁ = B[:,1]; b₂ = B[:,2]; b₃ = B[:,3]

function kdictGen(A)
	B = transpose(2*π*inv(A))
	kdict = Dict(
		"Γ" => B*[ 0  ;    0;   0],
		"Γ + iX₁" => B*[ im/2  ;    0;   0],
		"Γ + 10*iX₁" => B*[ 10*im/2  ;    0;   0],
		"A" => B*[ 1/2;  1/2; 1/2],
		"M" => B*[ 1/2;  1/2;   0],
		"Z" => B*[   0;    0; 1/2],
		"-Z" => B*[  0;    0;-1/2],
		"X₁" => B*[   1/2;    0; 0],
		"-X₁" => B*[  -1/2;    0;0],
		"X₂" => B*[   0;    1/2; 0],
		"-X₂" => B*[  0;  -1/2;0],
		"X₃" => B*[   0;    0; 1/2],
		"-X₃" => B*[  0;    0;-1/2],
		
		"-X₃" => B*[  0;    0;-1/2],
		"-X₃" => B*[  0;    0;-1/2],
		"-X₃" => B*[  0;    0;-1/2],
		"-X₃" => B*[  0;    0;-1/2],
		)
	return kdict
end

function βgen(p,runtype::String,β₀::Float64=0.2*eV, θ::Float64=360, startDWs::Float64=-5*nm)
	C = ħ/m₀
        β₀ = C^-1 * β₀
	if(runtype=="fmdotsP")
		function fmdotPβ(R::Vector{Float64})
			rad = 0.25; λ = 2*nm
			SLa₁ = p.A[:,1]; SLa₂ = p.A[:,2]
                        coeff = C^-1*β₀*exp(-((p.nz-1)*p.a₃[3] - R[3])/λ)*[1;0;0]
                        #println("Coeff = $(round.(C*coeff,sigdigits=3)), R = $R")
                        for i = 0:1
				for j = 0:1
					R₀ = i*SLa₁ + j*SLa₂
					if( (R[1]-R₀[1])^2 + (R[2]-R₀[2])^2 < (0.25*norm(SLa₂))^2)
						return 1*coeff
					end
				end
			end
			if( (R[1]-0.5*SLa₁[1])^2 + (R[2]-0.5*SLa₂[2])^2 < (0.25*norm(SLa₂))^2)
				return 1*coeff
			else
				return 0*coeff
			end
		end
                return fmdotPβ
	elseif(runtype=="multineeldws")
            function mndwβ(R::Vector{Float64})
                θ = 0
                error = 1
                # see: O'Handley on domain walls
                #a = 5*nm; l = 1*nm
                #a = 30*nm; l = 5*nm
                i = 0
                #decay = 1
                decay= exp(-((p.nz-1)*p.a₃[3] - R[3])/p.λ)
                while error > 10^-5
                    θ₋ = deepcopy(θ)
                    θ += 2*atan(exp(π*(p.startDWs - R[1] - p.DWspacing*i)/p.DWwidth))
                    i += 1
                    error = abs(θ-θ₋)
                end
                return β₀*[sin(θ);cos(θ);0]*decay
            end
            return mndwβ
	elseif(runtype=="multiblochdws")
            function mbdwβ(R::Vector{Float64})
                θ = 0
                error = 1
                # see: O'Handley on domain walls
                #a = 5*nm; l = 1*nm
                #a = 30*nm; l = 5*nm
                i = 0
                #decay = 1
                decay= exp(-((p.nz-1)*p.a₃[3] - R[3])/p.λ)
                while error > 10^-5
                    θ₋ = deepcopy(θ)
                    θ += 2*atan(exp(π*(p.startDWs - R[1] - p.DWspacing*i)/p.DWwidth))
                    i += 1
                    error = abs(θ-θ₋)
                end
                return β₀*[0;cos(θ);sin(θ)]*decay
            end
            return mbdwβ
	elseif(runtype=="fmdotsAP")
            function fmdotAPβ(R::Vector{Float64})
			rad = 0.25; λ = 2*nm
			SLa₁ = p.A[:,1]; SLa₂ = p.A[:,2]
                        coeff = C^-1*β₀*exp(-((p.nz-1)*p.a₃[3] - R[3])/λ)*[1;0;0]
                        #println("Coeff = $(round.(coeff,sigdigits=3)), R = $R")
                        for i = 0:1
				for j = 0:1
					R₀ = i*SLa₁ + j*SLa₂
					if( (R[1]-R₀[1])^2 + (R[2]-R₀[2])^2 < (0.25*norm(SLa₂))^2)
						return 1*coeff
					end
				end
			end
			if( (R[1]-0.5*SLa₁[1])^2 + (R[2]-0.5*SLa₂[2])^2 < (0.25*norm(SLa₂))^2)
				return -1*coeff
			else
				return 0*coeff
			end
		end
		return fmdotAPβ
	elseif(runtype=="neelwall")
		function dwnβ(R::Vector{Float64})
			α = 51 # arbitrary constant to smooth square wave
			λ = 2*nm # penetration depth of ferromagnetism into slab
                        xmag = cos(2*π*(θ/360)*(R⋅p.A[:,1]/(p.A[:,1]⋅p.A[:,1])))^α
                        ymag = (1-xmag^2)*sign(sin(2*π*(θ/360)*R⋅p.A[:,1]/(p.A[:,1]⋅p.A[:,1])))
			decay= C^-1*exp(-((p.nz-1)*p.a₃[3] - R[3])/λ)
			return β₀*[xmag;ymag;0]*decay
		end
		return dwnβ
	elseif(runtype=="simplejunction")
		function sjβ(R::Vector{Float64})
			α = 51 # arbitrary constant to smooth square wave
                        half = p.A[1,1]/2
			λ = 2*nm # penetration depth of ferromagnetism into slab
			decay= exp(-((p.nz-1)*p.a₃[3] - R[3])/λ)
                        if(R[1] > half*(1-0.00001))
                            return β₀*[0;cosd(θ);sind(θ)]*decay
                        else
                            return β₀*[0;1;0]*decay
                        end
		end
		return sjβ
	elseif(runtype=="blochdw")
		function dwbβ(R::Vector{Float64})
			α = 51 # arbitrary constant to smooth square wave
			λ = p.λ # penetration depth of ferromagnetism into slab
                        ymag = cos(2*π*(θ/360)*(R⋅p.A[:,1]/(p.A[:,1]⋅p.A[:,1])))
                        zmag = (1-ymag^2)*sign(sin(2*π*(θ/360)*R⋅p.A[:,1]/(p.A[:,1]⋅p.A[:,1])))
			decay=1
                        #decay= C^-1*exp(-((p.nz-1)*p.a₃[3] - R[3])/λ)
			return β₀*[0;ymag;zmag]*decay
		end
		return dwbβ
	elseif(runtype=="neellattice")
		function latnβ(R::Vector{Float64})
			α = 51 # arbitrary constant to smooth square wave
			λ = p.λ # penetration depth of ferromagnetism into slab
                        ymag = cos(2*π*(p.θ/360)*(R⋅p.A[:,1]/(p.A[:,1]⋅p.A[:,1])))
                        xmag = sin(2*π*(p.θ/360)*(R⋅p.A[:,1]/(p.A[:,1]⋅p.A[:,1])))
                        decay= exp(-((p.nz-1)*p.a₃[3] - R[3])/λ)
			return β₀*[xmag;ymag;0]*decay
		end
		return latnβ
	elseif(runtype=="blochlattice")
		function latbβ(R::Vector{Float64})
			α = 51 # arbitrary constant to smooth square wave
			λ = p.λ # penetration depth of ferromagnetism into slab
                        ymag = cos(2*π*(p.θ/360)*(R⋅p.A[:,1]/(p.A[:,1]⋅p.A[:,1])))
                        zmag = sin(2*π*(p.θ/360)*(R⋅p.A[:,1]/(p.A[:,1]⋅p.A[:,1])))
                        decay= exp(-((p.nz-1)*p.a₃[3] - R[3])/λ)
			return β₀*[0;ymag;zmag]*decay
		end
		return latbβ
	elseif(runtype=="fmthinfilm")
		function fmβ(R::Vector{Float64})
			λ = 5*nm
			decay= C^-1*exp(-((p.nz-1)*p.a₃[3] - R[3])/λ)
			return β₀*[0;1;0]*decay
		end
		return fmβ
	else
		function noβ(R::Vector{Float64})
			return [0.0;0.0;0.0]
		end
	end
end


params = (
	  t = t, t₁ = t₁, t₂ = t₂, t₃ = t₃, t₄ = t₄, t₅ = t₅, t₆ = t₆, t₇ = t₇, t₈ = t₈, t₉ = t₉,
	  vf = 10^6, η = 10^-4,
	  ε = ε, ε₁ = 2*eV, 
	  a₁ = a₁, a₂ = a₂, a₃ = a₃, A = A, a=a, b=a, c=c,
	  SLa₁ = a₁, SLa₂ = a₂, SLa₃ = a₃,
	  norb = 2, nsite = 1,
	  kdict = kdictGen(A), prune = [], μ_disorder = 0.0*eV
	  )

function pruneHoppingType(runtype::String="")
	println("runtype = $runtype")
	if(runtype ∈ ["nanopillars", "eggcarton", "afmthinfilm", "fmthinfilm", "fmdotsP", "fmdotsAP", "neelwall", "blochwall" ])
		return []
		#return ["z"]
	end
	if(runtype == "domainwallthin")
		return ["z","y"]
	end
	if(runtype == "device")
		return ["x","y","z"]
	end
	return []
end


function genLatSL(p,SL1::Vector{Int},SL2::Vector{Int},SL3::Vector{Int})
	SLa₁ = SL1[1]*p.a₁ + SL1[2]*p.a₂ + SL1[3]*p.a₃
	SLa₂ = SL2[1]*p.a₁ + SL2[2]*p.a₂ + SL2[3]*p.a₃
	SLa₃ = SL3[1]*p.a₁ + SL3[2]*p.a₂ + SL3[3]*p.a₃
	return SLa₁, SLa₂, SLa₃
end

function genSL(p,nx::Int,ny::Int,nz::Int,SL1::Vector{Int},SL2::Vector{Int},SL3::Vector{Int},runtype::String="",fieldtype="β")
	SLa₁, SLa₂, SLa₃ = genLatSL(p,SL1,SL2,SL3)
	newA = hcat(nx*p.SLa₁,ny*p.SLa₂,nz*p.SLa₃)
	if(nx*ny*nz*p.norb*p.nsite > 64)
		arpack = true
	else
		arpack = false
	end
        pruning = String["y","z"]
        SLparams = (
		SLa₁ = newA[:,1], SLa₂ = newA[:,2], SLa₃ = newA[:,3],
		A = newA, B = 2*π*transpose(inv(newA)),
		nx = nx, ny = ny, nz = nz, n=nx*ny*nz,
		kdict = kdictGen(newA),
		runtype=runtype,
		arpack=arpack,
		prune=pruning,
		klist = ["M","Γ","X₁","M","X₂","Γ", "X₃"],
		fieldtype=fieldtype,
		electrodeMaterial="metal",
		electrodeMagnetization = true,
		mixedDOS = false,
		nk = 0,
		E_samples = [i for i=-2.0:0.2:2.0],
		l_scattering = 0,
		δV = 0.01,
		T = 300,
		n_BLAS = 8,
		savedata = false,
		path = "../data-output/",
		startDWs = -5*nm, DWwidth = 9*nm, DWspacing = 15*nm, λ = 10^16*nm,
		β = 0.25*eV,
		η = 1*10^-4,
		ηD = 10^-4,
		θ=360.0
	)
	return merge(p,SLparams)
end

function generateParams()
	⊗(A,B) = kron(A,B)
	nx = 1; ny = 1; nz = 1;
    SL1 = [nx; 0; 0]; SL2 = [0; ny; 0]; SL3 = [0; 0; nz]

    runparams = (
        # fermi energy, σ of background disorder
        μ = 0.1*eV, μ_disorder = 0.0*eV, 

        # exchange splitting, name of field pattern, A or β (vector pot or exchange), finite-time broadenings
        runtype = "simplejunction", fieldtype = "β",

        # things to return 
        returnvals = ["transmission"],
        prune=[],
		deviceMaterial = "metal",
		deviceMagnetization = true
    )

    parms = merge(params,runparams)
    p = InBi.genSL(parms, nx, ny, nz, SL1, SL2, SL3, parms.runtype, parms.fieldtype) # generate SL params
    # supercell consisting of one unit cell
    pNNs = (n = p.n, nx = p.nx, ny = p.ny, nz = p.nz, nsite = p.nsite, norb = p.norb,
            t = 1, SLa₂ = p.SLa₂, a₁ = p.a₁, a₂ = p.a₂, a₃ = p.a₃, deviceMaterial = "metal", 
            ε₁ = p.ε₁, A = p.A, fieldtype = "A")
	pH = (μ = p.μ, μ_disorder = p.μ_disorder, n = p.n, norb = p.norb, nsite = p.nsite, 
			A = p.A, t = p.t, ε₁ = p.ε₁)
	pHopMat = (n = p.n, nsite = p.nsite, norb = p.norb, t = p.t, ε₁ = p.ε₁)
	A = βgen(p,p.runtype,p.β,p.θ,p.startDWs)

	println("Generated parameters in InBi")
	println(p)
    return p, pNNs, pH, pHopMat, A
end

end


