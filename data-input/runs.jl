# this is the top level definition of the run parameters

runparams = Dict(
	"path" => OUTPUT_DIR,
	"material_hamiltonian" => material_hamiltonians,
	"material_params" => Dict("t" => 1.0, "ε₀" => 1.0, "site_positions"=>site_positions),

	# define the routines to run. The three main ones are unitcell, transport, and supercell. 
	# Comment out a routine to not run it.

	# so, for runs looking at the electronic properties of just one unit cell
	"unitcell" => Dict("geometry" => devicegeometry, "material"=> "metal", "bands" => true, "bands_project" => [σ[1]], "save"=>[:bandstructure], "poisson" => false, "DOS" => false, "klist" => ["Γ","X₁","M","X₂","Γ","X₃"], "numInterpolations" => 64),
	
	# and for runs looking at the conductance of a large supercell
	"transport" => Dict("geometry" => devicegeometry, "ΔV" => 0.001, "μ" => 0.0*eV, "T" => 1, "η" => 1E-5*eV, "save" => [:transmission, :conductance], "electrodeMagnetization" => false, "deviceMagnetization" => false, "Gʳinv_method" => :RGF, "D_spin" => 0.000001*eV, "D_momentum" => 0.000005*eV, "kspace"=>false, "E_samples" => [E for E = -3.0:0.1:5.0], "electrodeMaterial" => "metal"),
	
	# and for runs where we want to slap a bunch of unit cells together and get the scattering-corrected electronic properties
	"supercell" => Dict("geometry" => devicegeometry, "bands_project" => [σ[1],σ[2]], "poisson"=>false, "μ" => 0.0*eV, "T" => 300, "η" => 1E-4*eV, "save" => [:unfoldedbands], "density_project" => [I(2),[σ[1],σ[2],σ[3]]], "Gʳinv_method" => :RGF, "D_dephasing" => 0.001*eV, "D_spin" => 0.0001*eV, "D_momentum" => 0.001*eV)
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
		merge!(runparams["unitcell"],geometry_params) # added this line for now
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
