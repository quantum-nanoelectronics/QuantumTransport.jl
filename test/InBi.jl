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

nx = ny = nz = 1; n = nx*ny*nz


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


params = (
	  t = t, t₁ = t₁, t₂ = t₂, t₃ = t₃, t₄ = t₄, t₅ = t₅, t₆ = t₆, t₇ = t₇, t₈ = t₈, t₉ = t₉,
	  vf = 10^6, η = 10^-3,
	  ε = ε, ε₁ = 2*eV, 
	  a₁ = a₁, a₂ = a₂, a₃ = a₃, A = A, a=a, b=a, c=c,
	  SLa₁ = a₁, SLa₂ = a₂, SLa₃ = a₃,
	  nx = nx, ny = ny, nz = nz, n = n, norb = 1, nsite = 1,
	  kdict = kdictGen(A), prune = [], μ_disorder = 0.025*eV
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

function genSL(p,nx::Int,ny::Int,nz::Int,SL1::Vector{Int},SL2::Vector{Int},SL3::Vector{Int},runtype::String="",fieldtype="A")
	SLa₁, SLa₂, SLa₃ = genLatSL(p,SL1,SL2,SL3)
	newA = hcat(nx*p.SLa₁,ny*p.SLa₂,nz*p.SLa₃)
	if(nx*ny*nz*p.norb*p.nsite > 64)
		arpack = true
	else
		arpack = false
	end
        pruning = String[]
        if(p.prune == [])
            pruning = pruneHoppingType(runtype)
        else
            pruning = p.prune
        end
        SLparams = (
		SLa₁ = newA[:,1], SLa₂ = newA[:,2], SLa₃ = newA[:,3],
		A = newA, B = 2*π*transpose(inv(newA)),
		nx = nx, ny = ny, nz = nz, n=nx*ny*nz,
		kdict = kdictGen(newA),
		runtype=runtype,
		arpack=arpack,
		prune=pruning,
		klist = ["M","Γ","X₁","M","X₂","Γ", "X₃"],
		fieldtype=fieldtype
	)
	return merge(p,SLparams)
end

function generateParams()
	⊗(A,B) = kron(A,B)
	nx = 1; ny = 1; nz = 1;
    SL1 = [nx; 0; 0]; SL2 = [0; ny; 0]; SL3 = [0; 0; nz]

    runparams = (
        # fermi energy, σ of background disorder
        μ = 0.0*eV, μ_disorder = 0.0*eV, 

        # exchange splitting, name of field pattern, A or β (vector pot or exchange), finite-time broadenings
        runtype = "blochlattice", fieldtype = "β",

        # things to return 
        returnvals = ["transmission"],
        prune=[]
    )

    parms = merge(params,runparams)
    p = InBi.genSL(parms, nx, ny, nz, SL1, SL2, SL3, parms.runtype, parms.fieldtype) # generate SL params
    # supercell consisting of one unit cell
    pNNs = (n = p.n, nx = p.nx, ny = p.ny, nz = p.nz, nsite = p.nsite, norb = p.norb,
            t = 1, SLa₂ = p.SLa₂, a₁ = p.a₁, a₂ = p.a₂, a₃ = p.a₃, deviceMaterial = "ins", 
            ε₁ = p.ε₁, A = p.A, fieldtype = "A")
	pH = (μ = p.μ, μ_disorder = p.μ_disorder, n = p.n, norb = p.norb, nsite = p.nsite, 
			A = p.A, t = p.t, ε₁ = p.ε₁)
	pHopMat = (n = p.n, nsite = p.nsite, norb = p.norb, t = p.t, ε₁ = p.ε₁)
    return pNNs, pH, pHopMat
end

end


