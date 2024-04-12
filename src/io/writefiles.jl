using CSV
using DataFrames
using Random
#using Tables
using Base.Filesystem: mktemp  # For creating a temporary file



function save_data_formatted(typeofdata::String, path::String, filename::String, axis_labels::Vector{String}, flip_axes::Bool=false, title::String="", data::Vector{Vector})
    if typeofdata == "ℝ→ℝ"
        x = data[1]
        y = data[2]
        # metadata is in the format ["type of plot", n_entries, n_codomain (dimensionality of codomain), title, flip_axes]
        metadata = [typeofdata, String(size(x)[1]), String(1), title, String(flip_axes)]
        df = DataFrame(
            Symbol(axis_labels[1]) = x,
            Symbol(axis_labels[2]) = y
        )
        save_csv(path, filename, df, metadata)
    elseif typeofdata == "ℝ³→R"
        x = data[1]
        y = data[2]
        z = data[3]
        C = data[4]
        # metadata is in the format ["type of plot", n_entries, n_codomain (dimensionality of codomain), title, flip_axes]
        metadata = [typeofdata, String(size(x)[1]), String(1), title, String(flip_axes)]
        df = DataFrame(
            Symbol(axis_labels[1]) = x,
            Symbol(axis_labels[2]) = y,
            Symbol(axis_labels[3]) = z,
            Symbol(axis_labels[4]) = C
        )
        save_csv(path, filename, df, metadata)
end

function save_csv(output_dir::String, filename::String, df::DataFrame, metadata::Vector{String})
    header = join(metadata, ",")

    # Create a temporary file to first write the DataFrame
    temp_csv_path, temp_file = mktemp()
    close(temp_file)  # Close the file because CSV.write opens it again

    # Define the final CSV file path with the optional output directory
    final_csv_path = joinpath(output_dir, filename)

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
