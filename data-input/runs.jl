### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 5c17c605-8bee-436a-8e9d-220e92c82e7e


# ╔═╡ ea227a19-a44e-479d-bfbe-650134d5d315


# ╔═╡ c9a1ae7a-faa3-11ee-1e5f-8b31657af949


# this is the top level definition of the run parameters
A_field(R::Vector{Float64}) = [0.0, 0.0, 0.0]

runparams = Dict(
	"path" => OUTPUT_DIR,
	"material_hamiltonian" => material_hamiltonians,
	"material_params" => Dict("t" => 1.0, "ε₀" => 1.0, "site_positions"=>site_positions),
    # Order of matrix_params: full matrix size, block size, phi, eta term, zeroThreshold term, σ₂, energy, 
	# (1000, 2, 0.2001, 1e-10, 1e-10, [0 -im; im 0], 3.0, matrixIndex)
	"matrix_params" => Dict("inv" => "RGF", "ϕ" => 0.0, "errorThreshold" => 1e-10),

	# A field
	"A_field" => A_field,

	# define the routines to run. The three main ones are unitcell, transport, and supercell. 
	# Comment out one of the dictionary entries below to not run the assocuated routine.

	# so, for runs looking at the electronic properties of just one unit cell
	"unitcell" => Dict("runtype" => "unitcell", "geometry" => devicegeometry, "material"=> "metal", "bands" => true, "bands_project" => [σ[1]], "save"=>[:bandstructure], "poisson" => false, "DOS" => false, "klist" => ["Γ","X₁","M","X₂","Γ","X₃"], "numInterpolations" => 64),
	
	# and for runs looking at the conductance of a large supercell
	"transport" => Dict("runtype" => "transport", "geometry" => devicegeometry, "ΔV" => 0.001, "μ" => 0.0*eV, "T" => 300, "η" => 1e-4, "save" => [:transmission, :conductance], "electrodeMagnetization" => false, "deviceMagnetization" => false, "Gʳinv_method" => :RGF, "D_spin" => 0.000001*eV, "D_momentum" => 0.000005*eV, "kspace"=>false, "E_samples" => [E for E = 0.0:1.0:10.0], "electrodeMaterial" => "metal"),

	# and for runs where we want to slap a bunch of unit cells together and get the scattering-corrected electronic properties
	"supercell" => Dict("runtype" => "supercell", "geometry" => devicegeometry, "bands_project" => [σ[1],σ[2]], "poisson"=>false, "μ" => 0.0*eV, "T" => 300, "η" => 1E-3*eV, "save" => [:unfoldedbands], "density_project" => [I(2),[σ[1],σ[2],σ[3]]], "Gʳinv_method" => :RGF, "D_dephasing" => 0.001*eV, "D_spin" => 0.0001*eV, "D_momentum" => 0.001*eV)
)

# all parameters of the simulation must be passed in through the runparams, and now we will ensure they are added
# or calculated in the beginning. 


function addCommonParams!(runparams)
	# this is a loop that will add keys to all routines
	for key in ["unitcell", "supercell", "transport"]
		if haskey(runparams,key)
			runparams[key]["material_hamiltonian"] = runparams["material_hamiltonian"]
			merge!(runparams[key], subspace_sizes)
			runparams[key]["path"] = runparams["path"]
			# Adding the A field function
			runparams[key]["A_field"] = runparams["A_field"]
		end
	end
end

function addTransportParams!(runparams)
	if haskey(runparams,"transport")
		params = runparams["transport"]
		# size of the H(k) hamiltonian for the entire device
		n_device = geometry_params["nx"]*geometry_params["ny"]*geometry_params["nz"]*subspace_sizes["nsite"]*subspace_sizes["norb"]*subspace_sizes["nspin"]
	
		# now we will add in the H(k) subspace, and then overwrite the nx, ny, nz with the device geometry
		merge!(params, runparams["matrix_params"], geometry_params, runparams["material_params"])

		params["n"] = n_device
		push!(params["prune"],"x")
		params["G"] = 2*π*pinv(params["A"])
		# runparams["transport"]["scattering"] = true
	end
end

function addUnitcellParams!(runparams)
	if haskey(runparams,"unitcell")
		# size of the H(k) hamiltonian for one unit cell
		n_unitcell = subspace_sizes["nsite"]*subspace_sizes["norb"]*subspace_sizes["nspin"]
		params = runparams["unitcell"]

		# added this line for now
		merge!(params, Dict("A"=>geometry_params["A"]), Dict("nx" => 1, "ny" => 1, "nz" => 1, "prune"=> []), runparams["material_params"]) 

		params["n"] = n_unitcell
	end
end

function addSupercellParams!(runparams)
	if haskey(runparams,"supercell")
		n_device = geometry_params["nx"]*geometry_params["ny"]*geometry_params["nz"]*subspace_sizes["nsite"]*subspace_sizes["norb"]*subspace_sizes["nspin"]
		params = runparams["supercell"]

		merge!(params, geometry_params, runparams["material_params"])
		
		params["n"] = n_device
	end
end

addCommonParams!(runparams)
addTransportParams!(runparams)
addUnitcellParams!(runparams)
addSupercellParams!(runparams)

# ╔═╡ Cell order:
# ╠═5c17c605-8bee-436a-8e9d-220e92c82e7e
# ╠═ea227a19-a44e-479d-bfbe-650134d5d315
# ╠═c9a1ae7a-faa3-11ee-1e5f-8b31657af949
