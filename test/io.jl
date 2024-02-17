using QuantumTransport
using Test
using LinearAlgebra
using DataFrames

# test the files in the io folder, but first generate positions


# contains constants used in make-cnt-data and other data visualization
const nm = 1E-9
const δ₀ = 0.142*nm
const a = 0.246*nm

const unicodeMeta = ["ħ", "ψ" , "σₓ", "ϕ"]


function Rx(θ::Float64) # returns a rotation matrix that rotates in the y-z plane, or around x
	A = [1 0 0; 0 cos(θ) -sin(θ); 0 sin(θ) cos(θ)]
	return A
end

function make_metallic_CNT_positions(nx::Int, radius::Float64=3*nm)
	positions = Vector{Float64}[]
	# gotta do some geometry
	circumferential_spacing = δ₀*3
	# number of graphene unit cells that make up CNT, circumferentially
	n_circ = Int(round(2*π*radius/circumferential_spacing))
	true_radius = n_circ*circumferential_spacing/(2*π) # maybe this is of use?

	# make the positions of the A and B sites in graphene
	posA = [0;-radius;δ₀/2]; posB = [0;-radius; -δ₀/2]
	# now we will drag them around and spin them around the nanotube to tile over space
	for ix = 0:(nx-1) # indexes the position along the nanotube
		for iθ = 0:(n_circ-1) # indexes position around the nanotube
			# okay. We will plop down 4 atoms per instance here
			θ = 2*π*iθ/n_circ # first two atoms
			Δx = ix*a
			atom1 = [Δx; 0; 0] + Rx(θ)*posA
			atom2 = [Δx; 0; 0] + Rx(θ)*posB
			# next two
			θ = 2*π*(iθ+0.5)/n_circ
			Δx = (ix+1/2)*a
			atom3 = [Δx; 0; 0] + Rx(θ)*posA
			atom4 = [Δx; 0; 0] + Rx(θ)*posB
			append!(positions,[atom1,atom2,atom3,atom4])
		end
	end
	return positions
end

CNT_positions = make_metallic_CNT_positions(20,1*nm)

const SAMPLE_POSITIONS = [
    [1.0, 0.0, 0.0],
    [2.0, 0.0, 0.0],
    [3.0, 0.0, 0.0],
    [4.0, 0.0, 0.0],
    [5.0, 0.0, 0.0],
    [6.0, 0.0, 0.0],
    [7.0, 0.0, 0.0],
    [8.0, 0.0, 0.0],
    [9.0, 0.0, 0.0],
    [10.0, 0.0, 0.0]
]

positions = CNT_positions
filename = "scatterplot.csv"

# generate data
baseDir = abspath(joinpath(@__DIR__, ".."))
ioDir = joinpath(baseDir, "data-output")
df, meta = generate_csv(positions)

# saving data
@test save_csv(ioDir, filename, df, meta)
println("-Writing-")
println("DataFrame: ")
println(first(df, 5))
println("Metadata: ")
println(meta)

# reading data
vals = get_data(ioDir, filename)
@test !isnothing(vals[1])
println("-Reading-")
println("DataFrame: ")
println(first(vals[1], 5))
println("Metadata: ")
println(vals[2])




# repeat for another file with unicode metadata + headers
filename = "scatterplot-unicode.csv"
rename!(df, names(df) .=> Symbol.(unicodeMeta))

# saving data
@test save_csv(ioDir, filename, df, unicodeMeta)
println("-Writing-")
println("DataFrame: ")
println(first(df, 5))
println("Metadata: ")
println(meta)
@test isfile(joinpath(ioDir, filename))


# reading data
vals = get_data(ioDir, filename)
@test !isnothing(vals[1])
println("-Reading-")
println("DataFrame: ")
println(first(vals[1], 5))
println("Metadata: ")
println(vals[2])
@test isfile(joinpath(ioDir, filename))

