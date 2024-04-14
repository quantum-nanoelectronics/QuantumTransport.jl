### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 5c17c605-8bee-436a-8e9d-220e92c82e7e


# ╔═╡ ea227a19-a44e-479d-bfbe-650134d5d315


# ╔═╡ c9a1ae7a-faa3-11ee-1e5f-8b31657af949


# this is the top level definition of the run parameters


γ⁵ = 0 # ?? not sure. this can go in constants? # TODO Vivian
runparams = Dict(
	"path" => OUTPUT_DIR,
	"material_hamiltonian" => material_hamiltonians,
	"material_params" => Dict("t" => 1.0, "ε₀" => 1.0, "site_positions"=>site_positions),

	# define the routines to run. The three main ones are unitcell, transport, and supercell. 
	# Comment out a routine to not run it.

	# so, for runs looking at the electronic properties of just one unit cell
	"unitcell" => Dict("material"=> "metal", "bands" => true, "bands_project" => [σ[1],γ⁵], "save"=>[:bandstructure, :DOS], "poisson" => false, "DOS" => false),
	
	# and for runs looking at the voltage-dependent transport properties
	"transport" => Dict("geometry" => devicegeometry, "ΔV" => 0.01, "μ" => 0.1*eV, "T" => 300, "η" => 1E-4*eV, "save" => [:transmission, :conductance], "electrodeMagnetization" => false, "Gʳinv_method" => :RGF, "D_spin" => 0.01*eV, "D_momentum" => 0.5*eV, "kspace"=>false, "E_samples" => [E for E = 0.0:0.1:1.0], "electrodeMaterial" => "metal"),
	
	# and for runs where we want to slap a bunch of unit cells together and get the scattering-corrected electronic properties
	"supercell" => Dict("geometry" => devicegeometry, "bands_project" => [σ[1],σ[2]], "poisson"=>true, "μ" => 0.1*eV, "T" => 300, "η" => 1E-4*eV, "save" => [:unfoldedbands], "density_project" => [I(2),[σ[1],σ[2],σ[3]]], "Gʳinv_method" => :RGF, "D_dephasing" => 0.1*eV, "D_spin" => 0.01*eV, "D_momentum" => 0.5*eV)
)

# all parameters of the simulation must be passed in through the runparams, and now we will ensure they are added
# or calculated in the beginning. 
function add_more_params!(runparams)
	# size of the H(k) hamiltonian for one unit cell
	n_unitcell = subspace_sizes["nsite"]*subspace_sizes["norb"]*subspace_sizes["nspin"]
	# size of the H(k) hamiltonian for the entire device
	n_device = geometry_params["nx"]*geometry_params["ny"]*geometry_params["nz"]*subspace_sizes["nsite"]*subspace_sizes["norb"]*subspace_sizes["nspin"]


	# this is a loop that will add keys to all routines
	for key in ["unitcell", "supercell", "transport"]
		if !haskey(runparams,key)
			continue
		end
		runparams[key]["material_hamiltonian"] = runparams["material_hamiltonian"]
		merge!(runparams[key], subspace_sizes)
		runparams[key]["path"] = runparams["path"]
	end
	
	if haskey(runparams,"unitcell")
		runparams["unitcell"]["n"] = n_unitcell
		merge!(runparams["unitcell"], runparams["material_params"])
	end

	if haskey(runparams,"supercell")
		runparams["supercell"]["n"] = n_device
		merge!(runparams["supercell"],geometry_params)
		merge!(runparams["supercell"], runparams["material_params"])
	end

	if haskey(runparams,"transport")
		runparams["transport"]["n"] = n_device
		
		# now we will add in the H(k) subspace, and then overwrite the nx, ny, nz with the device geometry
		merge!(runparams["transport"],geometry_params)
		push!(runparams["transport"]["prune"],"x")
		merge!(runparams["transport"], runparams["material_params"])
		runparams["transport"]["G"] = 2*π*inv(runparams["transport"]["A"])
	end
end

add_more_params!(runparams)


# All of these have references in this code and are used.
# copied this from random, but the code needs these in order to complete run. these need to be moved to a place
# TODO Vivian
anotherparamdict = Dict(
    "a₁" => [1.0e-9, 0.0, 0.0],
    "a₂" => [0.0, 1.0e-9, 0.0],
    "a₃" => [0.0, 0.0, 1.0e-9],
	"deviceMaterial" => "metal",
	"deviceMagnetization" => true,
	"ε₁" => 2,
	"fieldtype" => "β",
	"A" => [1 0 0; 0 1 0; 0 0 1],
    "returnvals" => ["transmission"],
    "μ_disorder" => 0.0,
    "SLa₂" => [0.0, 1.0e-9, 0.0],
    "δV" => 0.01,
    "l_scattering" => 0.0,
    "vf" => 1000000,
    "n_BLAS" => 8,
    "γ" => 0,
)
merge!(runparams["transport"], anotherparamdict)

# ╔═╡ Cell order:
# ╠═5c17c605-8bee-436a-8e9d-220e92c82e7e
# ╠═ea227a19-a44e-479d-bfbe-650134d5d315
# ╠═c9a1ae7a-faa3-11ee-1e5f-8b31657af949
