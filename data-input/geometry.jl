function devicegeometry(R::Vector{Float64})
	x = R[1]; y = R[2]; z = R[3];
	if (x<10*nm || x>50*nm)
		return "ins"
	else
		return "metal"
	end
end

geometry_params = Dict("A" => 2.866*nm*I(3), # this is the matrix of unit cell lattice vectors
		   "nx" => 50, "ny" => 5, "nz" => 1, # number of times to tile this cell over space in each direction to make whole device
		   "prune" => ["x","y","z"], # by default, system will set up periodic boundary conditions. This clips those hoppings
		   )
