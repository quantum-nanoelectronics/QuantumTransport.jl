### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 5c17c605-8bee-436a-8e9d-220e92c82e7e


# ╔═╡ ea227a19-a44e-479d-bfbe-650134d5d315


# ╔═╡ c9a1ae7a-faa3-11ee-1e5f-8b31657af949



function devicegeometry(R::Vector{Float64})
	x = R[1]; y = R[2]; z = R[3];
	return "metal"
	#=if (x<10*nm || x>50*nm)
		return "ins"
	else
		return "metal"
	end=#
end

geometry_params = Dict("A" => 2.866*nm*I(3), # this is the matrix of unit cell lattice vectors
		   "nx" => 2, "ny" => 1, "nz" => 1, # number of times to tile this cell over space in each direction to make whole device
		 	"prune" =>[]  
		   #"prune" => ["z","y"], # by default, system will set up periodic boundary conditions. This clips those hoppings
		   )

# ╔═╡ Cell order:
# ╠═5c17c605-8bee-436a-8e9d-220e92c82e7e
# ╠═ea227a19-a44e-479d-bfbe-650134d5d315
# ╠═c9a1ae7a-faa3-11ee-1e5f-8b31657af949
