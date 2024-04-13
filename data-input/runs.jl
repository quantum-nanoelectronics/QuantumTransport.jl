include("geometry.jl")

# The include statement for materials.jl was moved to common
include("materials.jl")

# this is the top level definition of the run parameters
runparams = Dict(
	"path" => "./",
	#"material_hamiltonian" => material_hamiltonians,
	"material_params" => Dict("t" => 1.0, "ε₀" => 1.0, "site_positions"=>site_positions),
	# so, for runs looking at the electronic properties of just one unit cell
	#  "unitcell" => Dict("material"=> "metal", "bands" => true, "bands_project" => [σ[1],γ⁵], "save"=>[:bandstructure, :DOS], "poisson" => false, "DOS" => false),
	# and for runs looking at the voltage-dependent transport properties
	"transport" => Dict("geometry" => devicegeometry, "ΔV" => 0.01, "μ" => 0.1*eV, "T" => 300, "η" => 1E-4*eV, "save" => [:transmission, :conductance], "electrodeMagnetization" => false,
	"Gʳinv_method" => :RGF, "D_spin" => 0.01*eV, "D_momentum" => 0.5*eV, "kspace"=>false, "E_samples" => [E for E = 0.0:0.1:2.0], "electrodeMaterial" => "metal"),
	# and for runs where we want to slap a bunch of unit cells together and get the scattering-corrected electronic properties
	#  "supercell" => Dict("geometry" => devicegeometry, "bands_project" => [σ[1],σ[2]], "poisson"=>true, "μ" => 0.1*eV, "T" => 300, "η" => 1E-4*eV, "save" => [:unfoldedbands], "density_project" => [I(2),[σ[1],σ[2],σ[3]]],
	#   "Gʳinv_method" => :RGF, "D_dephasing" => 0.1*eV, "D_spin" => 0.01*eV, "D_momentum" => 0.5*eV)
	)

anotherparamdict = Dict(
	# added this for GenBZ
	"G" => Diagonal([π*1e9, π*1e9, π*1e9]),
	# added this for path
	"path" => OUTPUT_DIR, 

	# TODO Vivian
	# copied this from random, but the code needs these in order to complete run. these need to be moved to a place
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

# all parameters of the simulation must be passed in through the runparams, and now we will ensure they are added
# or calculated in the beginning. 
function add_more_params!(runparams)
	# size of the H(k) hamiltonian for one unit cell
 	# n_unitcell = subspace_sizes["nsite"]*subspace_sizes["norb"]*subspace_sizes["nspin"]

	# runparams["unitcell"] = merge(subspace_sizes,runparams["unitcell"])
	if(haskey(runparams,"transport"))
		runparams["transport"]["path"] = runparams["path"]
		# runparams["unitcell"]["n"] = n_unitcell

		# now we will add in the H(k) subspace, and then overwrite the nx, ny, nz with the device geometry
		runparams["transport"] = merge(subspace_sizes,runparams["transport"])
		# runparams["transport"] = merge(subspace_sizes,runparams["supercell"])
		# overwriting
		merge!(runparams["transport"],geometry_params)
		# merge!(runparams["supercell"],geometry_params)

		n_device = geometry_params["nx"]*geometry_params["ny"]*geometry_params["nz"]*subspace_sizes["nsite"]*subspace_sizes["norb"]*subspace_sizes["nspin"]
		runparams["transport"]["n"] = n_device
		# runparams["supercell"]["n"] = n_device
		push!(runparams["transport"]["prune"],"x")

		# merge!(runparams["unitcell"], runparams["material_params"])
		merge!(runparams["transport"], runparams["material_params"])
	end
	# merge!(runparams["supercell"], runparams["material_params"])
	# runparams["transport"]["G"] = 2*π*inv(runparams["transport"]["A"])

end

runparams = add_more_params!(runparams)

merge!(runparams, anotherparamdict)
