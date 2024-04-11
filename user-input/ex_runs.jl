include("./ex_geometry.jl")
include("./ex_materials.jl")

runparams = Dict(
	     "path" => "./",
		 "material_hamiltonian" => material_hamiltonians,
	     # so, for runs looking at the electronic properties of just one unit cell
	     "unitcell" => Dict("bands" => true, "bands_project" => [σ[1],γ⁵], "save"=>[:bandstructure, :DOS], "poisson" => false, "DOS" => false),
	     # and for runs looking at the voltage-dependent transport properties
	     "transport" => Dict("nx"=>1, "ny" => 1, "nz"=>1, "geometry" => devicegeometry, "ΔV" => 0.01, "μ" => 0.1*eV, "T" => 300, "η" => 1E-4*eV, "save" => [:transmission, :conductance],
			  "Gʳinv_method" => :RGF, "D_spin" => 0.01*eV, "D_momentum" => 0.5*eV, "kspace"=>false),
	     # and for runs where we want to slap a bunch of unit cells together and get the scattering-corrected electronic properties
	     "supercell" => Dict("nx"=>1, "ny" => 1, "nz"=>1, "geometry" => devicegeometry, "bands_project" = [σ[1],σ[2]], "poisson"=>true, "μ" => 0.1*eV, "T" => 300, "η" => 1E-4*eV, "save" => [:unfoldedbands], "density_project" => [I(2),[σ[1],σ[2],σ[3]]],
			  "Gʳinv_method" => :RGF, "D_dephasing" => 0.1*eV, "D_spin" => 0.01*eV, "D_momentum" => 0.5*eV)
	     )

runparams = add_more_params(runparams)

calculate(runparams)
