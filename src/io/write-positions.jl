using CSV
using DataFrames
using Random

include("../data-vis/make-cnt-data.jl")


# This function is currently hard-coded to generate a df with random electron density and for testing purposes.
function generate_csv(data::Vector{Vector{Float64}})
    num_positions = length(data)
    randomElectronDensity = rand(1.0:100.0, num_positions)

    # Preallocate and directly construct DataFrame with all columns
    df = DataFrame(
        Matrix(hcat(data...))',  # Transpose the position matrix to match DataFrame format
        [:X, :Y, :Z]
    )
    df[!, :D] = randomElectronDensity  # Add the density column to the DataFrame

    metadata = ["metadata1", "metadata2"]

    return df, metadata
end

function insert_header(file_path::String, header::String)
    # Read the content of the file into an array of strings
    lines = readlines(file_path)

    # Insert a string at the second position
    if length(lines) >= 1
        insert!(lines, 2, header)
    else
        push!(lines, header)
    end

    # Write the updated content back to the file
    open(file_path, "w") do file
        for line in lines
            println(file, line)
        end
    end
end

function save_csv(output_dir::String, df::DataFrame, metadata::Vector{String})

    header = join(metadata, ",")
    # Define the CSV file path with the optional output directory
    csv_file_path = joinpath(output_dir, "scatterplot.csv")

    try
        CSV.write(csv_file_path, df)
        insert_header(csv_file_path, header)

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
