using CSV
using DataFrames
using Random
using Base.Filesystem: mktemp  # For creating a temporary file

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

function save_csv(output_dir::String, df::DataFrame, metadata::Vector{String})
    header = join(metadata, ",")

    # Create a temporary file to first write the DataFrame
    temp_csv_path, temp_file = mktemp()
    close(temp_file)  # Close the file because CSV.write opens it again

    # Define the final CSV file path with the optional output directory
    final_csv_path = joinpath(output_dir, "scatterplot.csv")

    try
        # Write the DataFrame to the temporary file
        CSV.write(temp_csv_path, df)

        # Open the final CSV file for writing
        open(final_csv_path, "w") do file
            # Write the header first
            println(file, header)
            # Then write the content of the temporary CSV file
            write(file, read(temp_csv_path))
        end

        println("CSV file '$final_csv_path' has been created.")
    catch e
        println("An IO error probably occurred: ", e)
        return false
    finally
        # Clean up: remove the temporary file
        rm(temp_csv_path; force = true)
    end
    return true
end
