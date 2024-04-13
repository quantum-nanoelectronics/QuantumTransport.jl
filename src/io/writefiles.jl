using CSV
using DataFrames
using Random
#using Tables
using Base.Filesystem: mktemp  # For creating a temporary file



function save_data_formatted(typeofdata::String, path::String, filename::String, axis_labels::Vector{String}, data::Vector{V}; row_input::Bool=false, flip_axes::Bool=false, title::String="") where V <: AbstractVector
    if row_input
        preproc = reduce(vcat,transpose.(positions))
        data = [preproc[:,i] for i = 1:size(preproc)[2]]
    end
    if typeofdata == "ℝ→ℝ"
        # metadata is in the format ["type of plot", n_entries, n_codomain (dimensionality of codomain), title, flip_axes]
        metadata = [typeofdata, string(size(data[1])[1]), string(1), title, string(flip_axes)]
        df = DataFrame(x = data[1], y = data[2])
    elseif typeofdata == "ℝ³→ℝ"
        x = data[1]
        y = data[2]
        z = data[3]
        C = data[4]
        # metadata is in the format ["type of plot", n_entries, n_codomain (dimensionality of codomain), title, flip_axes]
        metadata = [typeofdata, string(size(x)[1]), string(1), title, string(flip_axes)]
        df = DataFrame(x = x, y = y, z=z, C=C)
    end
    if(@isdefined df)
        rename!(df, Symbol.(axis_labels))
        save_csv(path, filename, df, metadata)
    else
        println("Saving failed")
    end
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
