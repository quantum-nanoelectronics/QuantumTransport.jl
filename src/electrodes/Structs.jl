# This file includes all the structs used for the module

# TODO Must specify the type of it here
struct AtomNN
	A1 # atom 1
	A2 # atom 2 
	δ # radius from A1 to A2
	N # unit cell of second atom, [Nx, Ny, Nz]⋅[a₁, a₂, a₃]
end


mutable struct Electrode
	xrange::Vector{Int}
	yrange::Vector{Int}
	zrange::Vector{Int}
	n::Int
	connectfrom::String # direction connecting fThreads.@threadsrom (i.e, "-x","+x",etc)
	type::String # material that the electrode is made of
	A::Function # magnitude of the exchange field
end

# # Also in Materials.jl
mutable struct Hopping
	a::Int # orbital/site index 1
	b::Int # orbital/site index 2 with PBC
	ia # index vector of site A
	ib # index vector of site B without PBC
	ra # location of atom a
	rb # location of atom b
	r # radius from a to b
	t  # hopping parameter affiliated with c†₂c₁ in spin basis. (i.e., t*I(2) or t*σ₊ might make sense)
	edge::Bool # does this hop off the edge of the superlattice?
	N # vector describing the [n₁;n₂;n₃]⋅[a₁;a₂;a₃] superlattice unit cell of site ib
	desc::String
end

