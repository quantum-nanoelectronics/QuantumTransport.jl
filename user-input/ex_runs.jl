include("./ex_geometry.jl")
include("./ex_materials.jl")


# this is the top level definition of the run parameters
runparams = Dict(
	     "path" => "./",
		 "material_hamiltonian" => material_hamiltonians,
		 "material_params" => (t = 1.0, ε₀ = 1.0),
	     # so, for runs looking at the electronic properties of just one unit cell
	     "unitcell" => Dict("material"="metal", "bands" => true, "bands_project" => [σ[1],γ⁵], "save"=>[:bandstructure, :DOS], "poisson" => false, "DOS" => false),
	     # and for runs looking at the voltage-dependent transport properties
	     "transport" => Dict("geometry" => devicegeometry, "ΔV" => 0.01, "μ" => 0.1*eV, "T" => 300, "η" => 1E-4*eV, "save" => [:transmission, :conductance],
			  "Gʳinv_method" => :RGF, "D_spin" => 0.01*eV, "D_momentum" => 0.5*eV, "kspace"=>false),
	     # and for runs where we want to slap a bunch of unit cells together and get the scattering-corrected electronic properties
	     "supercell" => Dict("geometry" => devicegeometry, "bands_project" = [σ[1],σ[2]], "poisson"=>true, "μ" => 0.1*eV, "T" => 300, "η" => 1E-4*eV, "save" => [:unfoldedbands], "density_project" => [I(2),[σ[1],σ[2],σ[3]]],
			  "Gʳinv_method" => :RGF, "D_dephasing" => 0.1*eV, "D_spin" => 0.01*eV, "D_momentum" => 0.5*eV)
	     )

runparams = add_more_params!(runparams)

main(runparams)


# all parameters of the simulation must be passed in through the runparams, and now we will ensure they are added
# or calculated in the beginning. 
function add_more_params!(runparams)
	# size of the H(k) hamiltonian for one unit cell
 	n_unitcell = subspace_sizes["nsite"]*subspace_size["norb"]*subspace_sizes["nspin"]
	runparams["unitcell"] = merge(subspace_sizes,runparams["unitcell"])
	runparams["unitcell"]["n"] = n_unitcell

	# now we will add in the H(k) subspace, and then overwrite the nx, ny, nz with the device geometry
	runparams["transport"] = merge(subspace_sizes,runparams["transport"])
	runparams["supercell"] = merge(subspace_sizes,runparams["supercell"])
	
	# overwriting
	runparams["transport"] = merge(runparams["transport"],geometry_params)
	runparams["supercell"] = merge(runparams["supercell"],geometry_params)

	n_device = geometry_params["nx"]*geometry_params["ny"]*geometry_params["nz"]*subspace_sizes["nsite"]*subspace_sizes["norb"]*subspace_sizes["nspin"]
	runparams["transport"]["n"] = n_device
	runparams["supercell"]["n"] = n_device
	append!(runparams["transport"]["prune"],"x")
end
