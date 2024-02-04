using CSV
using DataFrames
using Random

include("../data-vis/make-cnt-data.jl")

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

function generate_csv(output_dir::String, positions::Vector{Vector{Float64}}=SAMPLE_POSITIONS)
    try
        num_positions = length(positions)
        randomElectronDensity = rand(1.0:100.0, num_positions)

        # Preallocate and directly construct DataFrame with all columns
        df = DataFrame(
            Matrix(hcat(positions...))',  # Transpose the position matrix to match DataFrame format
            [:X, :Y, :Z]
        )
        df[!, :D] = randomElectronDensity  # Add the density column

        # Define the CSV file path with the optional output directory
        csv_file_path = joinpath(output_dir, "scatterplot.csv")
        CSV.write(csv_file_path, df)

        println("CSV file '$csv_file_path' has been created.")
    catch e
        println("An IO error probably occurred: ", e)
        return false
    end
    return true
end


# positions = make_metallic_CNT_positions(20,1*nm)

# generate_csv("src/io/", positions)

# generate_csv("src/io/")
