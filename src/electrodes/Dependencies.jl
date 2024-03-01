include("Structs.jl")
include("NearestNeighbors.jl")
include("Utilities.jl")
include("Hamiltonians.jl")
include("SelfEnergies.jl")


# From Operators.jl 
σ = Vector{Matrix}(undef,3)
σ₀ = [
      1 0
      0 1
     ]

σ₁ = [
      0  1 
      1  0
     ] 
σ₂ = [
      0 -im 
      im  0
     ] 
σ₃ = [
      1  0
      0 -1
     ]
σ[1] = σ₁; σ[2] = σ₂; σ[3] = σ₃

⊗(A,B) = kron(A,B)
×(u,v) = cross(u,v)

# from Constants.jl
ħ = 1.05457E-34
h = ħ * 2*π
m₀ = 9.10938E-31
q = 1.60218E-19
ϵ₀ = 8.854E-12
metre = 1
au = 27.2
eV = 1.0
cm = 1E-2
μ₀ = 1.2566*10^-6 # H/m vacuum permeability
nm = 1E-9
Å = 1E-10
r₀ = 5.29E-11
Ry = m₀*q^4 / (8*h^2*ϵ₀^2)
μₑ = 9.28*10^-24 # electron magnetic moment in A*m^2
μB = 5.788838E-5 # bohr magneton in eV/T
kB = 8.617E-5 # boltzmann constant in eV/K


#From Materials.jl
function metalHopping(p::NamedTuple,NNs::Vector{Hopping},ia::Vector{Int})
	iorb = ia[5]
	t = (3*p.t - 1*eV)*(I(2))
	pushHopping!(NNs, t, ia, ia, p)
	for ax = 1:3
		for dir = [-1,1]
			# for weyl term in hamiltonian
			di = zeros(5); di[ax] = dir; ib = Int.(ia + di)
			#Ra = xyztor(p,ia); Rb = xyztor(p,ib); 
			#δ = Rb - Ra
			# implement H = +vf*𝐩⋅𝛔 = -vf𝑖ħ ∇ᵣ⋅σ on finite grid
                        t = -1*p.t*I(2)
			pushHopping!(NNs, t, ia, ib, p)
                        #t = (p.ϵ₁ + 2*p.t)*(I(2))
                        #pushHopping!(NNs, t, ia, ia, p)
		end
	end
end

function insHopping(p::NamedTuple,NNs::Vector{Hopping},ia::Vector{Int})
	iorb = ia[5]
	t = (p.ϵ₁ + 3*p.t)*(I(2))
	pushHopping!(NNs, t, ia, ia, p)
	for ax = 1:3
		for dir = [-1,1]
			# for weyl term in hamiltonian
			di = zeros(5); di[ax] = dir; ib = Int.(ia + di)
                        t = -1/2*p.t*I(2)
			pushHopping!(NNs, t, ia, ib, p)
                        #t = (p.ϵ₁ + 2*p.t)*(I(2))
                        #pushHopping!(NNs, t, ia, ia, p)
		end
	end
end

function nextsite(iorb::Int)
    return (1-2*iorb)
end


function pushHopping!(NNs::Vector, t, ia::Vector{Int}, ib::Vector{Int}, p) 
	a = xyztoi(p,ia); b = xyztoi(p,ib);
	ra = xyztor(p,ia); rb = xyztor(p,ib); r = rb - ra;
	# for hopping term
	NN = deepcopy(Hopping(a,b,ia,ib,ra,rb,r,t, false, [0;0;0],""))
	push!(NNs,NN)
end

function weyl3Hopping(p::NamedTuple,NNs::Vector{Hopping},ia::Vector{Int})
        iorb = ia[5]
        ib = ia; 
	# nearest neighbors
        for ax = 1:3
		for dir = [-1,1]
			# for weyl term in hamiltonian
			di = zeros(5); di[ax] = dir; 
                        #di[5] = nextsite(iorb); 
                        ib = Int.(ia + di)
                        t = zeros(ComplexF64,2,2)
                        if(ax == 3)
                            t = p.t*σ[3]
                        elseif(ax==1)
                            t == p.t*(-σ[3] .+ dir*2*im*σ[1])
                        elseif(ax==2)
                            t == p.t*(-σ[3] .+ dir*2*im*σ[2])
                        else
                            println("Something broken! No 4th axis")
                        end
			pushHopping!(NNs, t, ia, ib, p)
		end
	end
        # next nearest neighbors
        for dx = [-1,1]
            for dy = [-1,1]
                di = zeros(5); di[1] = dx; di[2] = dy;
                ib = Int.(ia+di)
                t = 1.5*im*p.t*(-dx*σ[1] .+ -dy*σ[2]) # implements the next-nearest neigbor term
                pushHopping!(NNs, t, ia, ib, p)
            end
        end
        # next-to-next nearest neighbors
        for ax = 1:2
            for dir = [-1,1]
                di = zeros(5); di[ax] = 2*dir
                ib = Int.(ia+di)
                t = dir*p.t*0.5*im*σ[ax]
                pushHopping!(NNs, t, ia, ib, p)
            end
        end
end


function weyl2Hopping(p::NamedTuple,NNs::Vector{Hopping},ia::Vector{Int})
        iorb = ia[5]
        ib = ia; 
	# nearest neighbors
        for ax = 1:3
		for dir = [-1,1]
			# for weyl term in hamiltonian
			di = zeros(5); di[ax] = dir; 
                        #di[5] = nextsite(iorb); 
                        ib = Int.(ia + di)
                        t = zeros(ComplexF64,2,2)
                        if(ax == 1)
                            t = p.t*(σ[1].-σ[3])
                        elseif(ax==2)
                            t = -p.t*(σ[1].+σ[3])
                        elseif(ax==3)
                            t = p.t*σ[3]
                        else
                            println("Something broken! No 4th axis")
                        end
                        #t = (-im/2)*dir*p.t*σ[ax]
			pushHopping!(NNs, t, ia, ib, p)
		end
	end
        for dx = [-1,1]
            for dy = [-1,1]
                di = zeros(5); di[1] = dx; di[2] = dy;
                ib = Int.(ia+di)
                t = -dx*dy*0.5*p.t*σ[2] # implements the next-nearest neigbor term
                pushHopping!(NNs, t, ia, ib, p)
            end
        end
end

function chern2DHopping(p::NamedTuple,NNs::Vector{Hopping},ia::Vector{Int})
        iorb = ia[5]
        ib = ia; 
        #ib[5] += nextsite(iorb)
        #t = 3*p.t*(I(2))
	#pushHopping!(NNs, t, ia, ib , p)
        # for H₂ = τ₃⊗σ₀
        t = nextsite(iorb)*(2*p.m2+p.γ)*(I(2))
	pushHopping!(NNs, t, ia, ia, p)
	for ax = 1:3
		for dir = [-1,1]
			# for weyl term in hamiltonian
			di = zeros(5); di[ax] = dir; 
                        # for Hweyl = τ₁⊗k⋅σ
                        di[5] = nextsite(iorb); 
                        ib = Int.(ia + di)
                        t = (im/2)*dir*p.t*σ[ax]
			#Ra = xyztor(p,ia); Rb = xyztor(p,ib); 
			#δ = Rb - Ra
			# implement H = +vf*𝐩⋅𝛔 = -vf𝑖ħ ∇ᵣ⋅σ on finite grid
                        #t = nextsite(iorb)*(-im/2)*dir*p.t*σ[ax]
			pushHopping!(NNs, t, ia, ib, p)
			# for normal hopping term in hamiltonian
			ib[5] = iorb; 
			
			t = -(1/2)*nextsite(iorb)*p.m2*(I(2))
			pushHopping!(NNs, t, ia, ib, p)
		end
	end
end


function weylHopping(p::NamedTuple,NNs::Vector{Hopping},ia::Vector{Int})
        iorb = ia[5]
        ib = ia; 
        #ib[5] += nextsite(iorb)
        #t = 3*p.t*(I(2))
	#pushHopping!(NNs, t, ia, ib , p)
        # for H₂ = τ₃⊗σ₀
        t = 3*nextsite(iorb)*p.t*(I(2))
	pushHopping!(NNs, t, ia, ia, p)
	for ax = 1:3
		for dir = [-1,1]
			# for weyl term in hamiltonian
			di = zeros(5); di[ax] = dir; 
                        # for Hweyl = τ₁⊗k⋅σ
                        di[5] = nextsite(iorb); 
                        ib = Int.(ia + di)
                        t = (-im/2)*dir*p.t*σ[ax]
			#Ra = xyztor(p,ia); Rb = xyztor(p,ib); 
			#δ = Rb - Ra
			# implement H = +vf*𝐩⋅𝛔 = -vf𝑖ħ ∇ᵣ⋅σ on finite grid
                        #t = nextsite(iorb)*(-im/2)*dir*p.t*σ[ax]
			pushHopping!(NNs, t, ia, ib, p)
			# for normal hopping term in hamiltonian
			ib[5] = iorb; 
			
			t = -(1/2)*nextsite(iorb)*p.t*(I(2))
			pushHopping!(NNs, t, ia, ib, p)
		end
	end
end


# from UsefulFunctions.jl

function genTetBZ(p::NamedTuple,nx::Int=0, ny::Int=100, nz::Int=100) # only works for cubic lattice
    # nx, ny, and nz specifically refer to # of points in IBZ
    kpoints = Vector{Float64}[]
    kindices = Vector{Int}[]
    kweights = Float64[]
    X1 = p.kdict["X₁"];
    X2 = p.kdict["X₂"];
    X3 = p.kdict["X₃"];
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
    return kpoints, kweights, kindices
end

hoppingDict = Dict{String,Function}(
                                    "mtjweyl"=>weylHopping,
                                    "2Dchern"=>chern2DHopping,
                                    "weyl"=>weylHopping,
                                    "weyl2"=>weyl2Hopping,
                                    "weyl3"=>weyl3Hopping,
                                    "wins"=>insHopping,
                                    "ins"=>insHopping,
                                    "insulator"=>insHopping,
                                    "metal"=>metalHopping
                                   )



