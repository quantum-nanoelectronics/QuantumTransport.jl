### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 4f25fe74-05f2-11ef-0956-959921a1d25e
using Pkg

# ╔═╡ e628282c-e5a5-4efe-bf59-7db4c30e0f49
Pkg.add(url="https://github.com/quantum-nanoelectronics/QuantumTransport.jl", rev="dev")

# ╔═╡ 9bb45983-9b5a-4215-8e91-72617e570815
Pkg.add("WGLMakie")

# ╔═╡ 55a58528-ff17-40cf-9f4d-9d6151225943
Pkg.add("Test")

# ╔═╡ 47d1595c-7963-4056-ba06-f325a88ed947
module testDriver

using QuantumTransport
using Test

function printDict(runparams, type)
    println("\033[1m$(type):\033[0m")
    for (key, value) in runparams
        println("  \033[1m$(key):\033[0m $(value)")
    end
end

function driverTest(params)
    main(params)
    return true
end


printDict(runparams["transport"], "Transport Parameters")
printDict(runparams["unitcell"], "Unitcell Parameters")
# printDict(runparams["supercell"], "Supercell Parameters")

@test driverTest(runparams)


end

# ╔═╡ 30f0ce34-020e-40a6-bd84-8ba3fca4a73d
using QuantumTransport

# ╔═╡ 667179fc-3898-4727-a797-6d0a9d7b7624
using GLMakie

# ╔═╡ bd129adc-74b1-48e2-afba-486dd326a199
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

    baseDir = abspath(joinpath(@__DIR__, ".."))
    ioDir = joinpath(baseDir, "data-output")
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


# ╔═╡ a9a5ff63-870b-4d58-af63-53f8790d8227
Pkg.rm("QuantumTransport")

# ╔═╡ 819b84fa-f5aa-4c71-b8cc-10bdb78b29ba
Pkg.update()

# ╔═╡ 01e5c475-0b94-4dc8-833e-c8df3adec628
# The following code is included in the test/testDriver.jl file

# ╔═╡ e92784a7-dc6d-4556-a87c-f74d81137b68
GLMakie.activate!()

# ╔═╡ 631d1bd6-be1e-49d6-9d7f-656a0d8e87ce
fig1 = call_function_based_on_header(OUTPUT_DIR, "transmission.csv")

# ╔═╡ f32522f2-5cc9-4cae-b182-57a71b01599f
fig2 = call_function_based_on_header(OUTPUT_DIR, "CNTpositions.csv")

# ╔═╡ a67082be-d618-4ac0-8bdb-beacb9aed4d9


# ╔═╡ bf400821-74bc-431f-918a-b5aea0eefd42


# ╔═╡ Cell order:
# ╠═4f25fe74-05f2-11ef-0956-959921a1d25e
# ╠═a9a5ff63-870b-4d58-af63-53f8790d8227
# ╠═819b84fa-f5aa-4c71-b8cc-10bdb78b29ba
# ╠═e628282c-e5a5-4efe-bf59-7db4c30e0f49
# ╠═9bb45983-9b5a-4215-8e91-72617e570815
# ╠═01e5c475-0b94-4dc8-833e-c8df3adec628
# ╠═47d1595c-7963-4056-ba06-f325a88ed947
# ╠═30f0ce34-020e-40a6-bd84-8ba3fca4a73d
# ╠═667179fc-3898-4727-a797-6d0a9d7b7624
# ╠═e92784a7-dc6d-4556-a87c-f74d81137b68
# ╠═631d1bd6-be1e-49d6-9d7f-656a0d8e87ce
# ╠═f32522f2-5cc9-4cae-b182-57a71b01599f
# ╠═a67082be-d618-4ac0-8bdb-beacb9aed4d9
# ╠═55a58528-ff17-40cf-9f4d-9d6151225943
# ╠═bd129adc-74b1-48e2-afba-486dd326a199
# ╠═bf400821-74bc-431f-918a-b5aea0eefd42