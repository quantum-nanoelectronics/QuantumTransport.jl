# This file is used to test the IO code

module TestIO

using QuantumTransport
using Test
using DataFrames
using Random
#using Makie #If this line is uncommented or if Makie is added to this package, github tests will fail


"""
	testIO(ioDir, filename, positions, meta, header)

This function performs the tests for this file.

# Arguments
- `ioDir`: The directory where the file will be saved.
- `filename`: The name of the file.
- `positions`: The positions of the elements.
- `meta`: The metadata associated with the file.
- `header`: The header information.

# Returns
- None
"""
function testIO(ioDir, filename, positions, meta, header)
    # generate data
    df, meta = generate_df(positions, meta, header)

    # saving data
    @test save_csv(ioDir, filename, df, meta)
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
end

function test_formatted_save(ioDir, filename)
    # generate data
    
    # saving data
    xvals = LinRange(0,1.0,10)
    yvals = rand(10)
    @test save_data_formatted("ℝ→ℝ", ioDir, filename, ["X (nm)", "Y (nm)"], [xvals,yvals], true)
    @test isfile(joinpath(ioDir, filename))


    # reading data
    #=vals = get_data(ioDir, filename)
    @test !isnothing(vals[1])
    println("-Reading-")
    println("DataFrame: ")
    println(first(vals[1], 5))
    println("Metadata: ")
    println(vals[2])
    @test isfile(joinpath(ioDir, filename))=#
end


"""
	make_metallic_CNT_positions(nx::Int, radius::Float64=3.0, a::Float64=0.246 * 1E-9, δ₀::Float64 = 0.142 * 1E-9)::Vector{Vector{Float64}}

Constructs the positions of metallic carbon nanotubes (CNTs) in a two-dimensional lattice.

# Arguments
- `nx::Int`: The number of CNTs in the x-direction.
- `radius::Float64`: The radius of the CNTs. Default is 3.0.
- `a::Float64`: The lattice constant. Default is 0.246 * 1E-9.
- `δ₀::Float64`: The displacement between neighboring CNTs. Default is 0.142 * 1E-9.

# Returns
A vector of vectors containing the positions of the metallic CNTs.

"""

function make_metallic_CNT_positions(nx::Int, radius::Float64=3.0 * 1E-9, a::Float64=0.246 * 1E-9, δ₀::Float64 = 0.142 * 1E-9)::Vector{Vector{Float64}}
	function Rx(θ::Float64) # returns a rotation matrix that rotates in the y-z plane, or around x
		A = [1 0 0; 0 cos(θ) -sin(θ); 0 sin(θ) cos(θ)]
		return A
	end

    positions = Vector{Float64}[]
    # gotta do some geometry
    circumferential_spacing = δ₀ * 3
    # number of graphene unit cells that make up CNT, circumferentially
    n_circ = Int(round(2 * π * radius / circumferential_spacing))
    true_radius = n_circ * circumferential_spacing / (2 * π) # maybe this is of use?

    # make the positions of the A and B sites in graphene
    posA = [0; -radius; δ₀ / 2]
    posB = [0; -radius; -δ₀ / 2]
    # now we will drag them around and spin them around the nanotube to tile over space
    for ix = 0:(nx-1) # indexes the position along the nanotube
        for iθ = 0:(n_circ-1) # indexes position around the nanotube
            # okay. We will plop down 4 atoms per instance here
            θ = 2 * π * iθ / n_circ # first two atoms
            Δx = ix * a
            atom1 = [Δx; 0; 0] + Rx(θ) * posA
            atom2 = [Δx; 0; 0] + Rx(θ) * posB
            # next two
            θ = 2 * π * (iθ + 0.5) / n_circ
            Δx = (ix + 1 / 2) * a
            atom3 = [Δx; 0; 0] + Rx(θ) * posA
            atom4 = [Δx; 0; 0] + Rx(θ) * posB
            append!(positions, [atom1, atom2, atom3, atom4])
        end
    end
    posmatrix = reduce(vcat,transpose.(positions))
    return [posmatrix[:,1],posmatrix[:,2],posmatrix[:,3]]
end


"""
	generate_df(data::Vector{Vector{Float64}}, metaData::Vector{String}, header)::Tuple{DataFrame,Vector{String}}

Generate a CSV file from the given data and metadata.

# Arguments
- `data::Vector{Vector{Float64}}`: The data to be written to the CSV file.
- `metaData::Vector{String}`: The metadata associated with the data.
- `header`: Header for the CSV file. If empty, the header will remain unchanged.

# Returns
- `Tuple{DataFrame,Vector{String}}`: A tuple containing the DataFrame representing the CSV file and the header.

"""
function generate_df(data::Vector{Vector{Float64}}, metaData::Vector{String}, header)::Tuple{DataFrame,Vector{String}}
    num_positions = length(data)
    randomElectronDensity = rand(1.0:100.0, num_positions)

    # Preallocate and directly construct DataFrame with all columns
    df = DataFrame(
        Matrix(hcat(data...))',  # Transpose the position matrix to match DataFrame format
        [:X, :Y, :Z]
    )
    df[!, :D] = randomElectronDensity  # Add the density column to the DataFrame

    if isempty(header)
        return df, metaData
    end

    # Rename the header of the DataFrame
    rename!(df, names(df) .=> Symbol.(header))
    return df, metaData
end

sample_positions = [
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

# set data and call test functions
function runIOTests()
    println("\033[1mGenerating CNT positions\033[0m")
    N = 20
    CNT_positions = make_metallic_CNT_positions(N)

    println("\033[1mGenerated CNT positions\033[0m")

    #=positions = CNT_positions
    filename1 = "scatterplot.csv"
    filename2 = "scatterplot-unicode.csv"
    metadata = ["metadata1", "metadata2"]
    unicodeMeta = ["ħ", "ψ" , "σₓ", "ϕ"]
    header = copy(unicodeMeta)
    testIO(ioDir, filename1, positions, metadata, [])
    testIO(ioDir, filename2, positions, unicodeMeta, header)=#
    xvals = LinRange(0,1,12)
    yvals = rand(12)

    @test save_data_formatted("ℝ→ℝ", OUTPUT_DIR, "testvals.csv", ["X (nm)", "Y (nm)"], [xvals,yvals])
    @test isfile(joinpath(OUTPUT_DIR, "testvals.csv"))


    @test save_data_formatted("ℝ³→ℝ", OUTPUT_DIR, "CNTpositions.csv", 
    ["X (nm)", "Y (nm)", "Z (nm)", "n (# electrons)"], [CNT_positions[1],CNT_positions[2],CNT_positions[3], rand(size(CNT_positions[1])[1])])
    @test isfile(joinpath(OUTPUT_DIR, "CNTpositions.csv"))
    #test_formatted_save(ioDir, "testvals.csv", [xvals, yvals])
    #test_formatted_save(ioDir, "CNT_positions.csv", [CNT_positions[1],CNT_positions[2],CNT_positions[3],])

end

runIOTests()

end # module MyModule
