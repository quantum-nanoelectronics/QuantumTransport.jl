### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 5c17c605-8bee-436a-8e9d-220e92c82e7e


# ╔═╡ ea227a19-a44e-479d-bfbe-650134d5d315


# ╔═╡ c9a1ae7a-faa3-11ee-1e5f-8b31657af949


# this is the top level definition of the run parameters

runparams = Dict(
	"path" => OUTPUT_DIR,
	"material_hamiltonian" => material_hamiltonians,
	"material_params" => Dict("t" => 1.0, "ε₀" => 1.0, "site_positions"=>site_positions),
    # Order of matrix_params: full matrix size, block size, phi, eta term, zeroThreshold term, σ₂, energy, 
	# (1000, 2, 0.2001, 1e-10, 1e-10, [0 -im; im 0], 3.0, matrixIndex)
	"matrix_params" => Dict("inv" => "RGF", "ϕ" => 0.2001, "errorThreshold" => 1e-10),

	# define the routines to run. The three main ones are unitcell, transport, and supercell. 
	# Comment out a routine to not run it.

	# so, for runs looking at the electronic properties of just one unit cell
	"unitcell" => Dict("geometry" => devicegeometry, "material"=> "metal", "bands" => true, "bands_project" => [σ[1]], "save"=>[:bandstructure], "poisson" => false, "DOS" => false, "klist" => ["Γ","X₁","M","X₂","Γ","X₃"], "numInterpolations" => 64),
	
	# and for runs looking at the conductance of a large supercell
	"transport" => Dict("geometry" => devicegeometry, "ΔV" => 0.001, "μ" => 0.0*eV, "T" => 300, "η" => 1E-2*eV, "save" => [:transmission, :conductance], "electrodeMagnetization" => false, "deviceMagnetization" => false, "Gʳinv_method" => :RGF, "D_spin" => 0.000001*eV, "D_momentum" => 0.000005*eV, "kspace"=>false, "E_samples" => [E for E = -0.2:0.025:1.4], "electrodeMaterial" => "metal"),
	
	# and for runs where we want to slap a bunch of unit cells together and get the scattering-corrected electronic properties
	"supercell" => Dict("geometry" => devicegeometry, "bands_project" => [σ[1],σ[2]], "poisson"=>false, "μ" => 0.0*eV, "T" => 300, "η" => 1E-3*eV, "save" => [:unfoldedbands], "density_project" => [I(2),[σ[1],σ[2],σ[3]]], "Gʳinv_method" => :RGF, "D_dephasing" => 0.001*eV, "D_spin" => 0.0001*eV, "D_momentum" => 0.001*eV)
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
		merge!(runparams["unitcell"], Dict("A"=>geometry_params["A"])) # added this line for now
		merge!(runparams["unitcell"], Dict("nx" => 1, "ny" => 1, "nz" => 1, "prune"=> [])) # added this line for now
		merge!(runparams["unitcell"], runparams["material_params"])
	end

	if haskey(runparams,"supercell")
		runparams["supercell"]["n"] = n_device
		merge!(runparams["supercell"],geometry_params)
		merge!(runparams["supercell"], runparams["material_params"])
	end

	if haskey(runparams,"transport")
		runparams["transport"]["n"] = n_device

		# runparams["transport"]["scattering"] = true
		merge!(runparams["transport"], runparams["matrix_params"])
		
		# now we will add in the H(k) subspace, and then overwrite the nx, ny, nz with the device geometry
		merge!(runparams["transport"],geometry_params)
		push!(runparams["transport"]["prune"],"x")
		merge!(runparams["transport"], runparams["material_params"])
		runparams["transport"]["G"] = 2*π*pinv(runparams["transport"]["A"])
	end
end

add_more_params!(runparams)

# ╔═╡ Cell order:
# ╠═5c17c605-8bee-436a-8e9d-220e92c82e7e
# ╠═ea227a19-a44e-479d-bfbe-650134d5d315
# ╠═c9a1ae7a-faa3-11ee-1e5f-8b31657af949
