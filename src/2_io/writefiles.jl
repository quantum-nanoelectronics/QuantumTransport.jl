using CSV
using DataFrames
using Random
using Serialization
using Base64
using Base.Filesystem: mktemp  # For creating a temporary file

"""
    save_data(typeofdata::Symbol, path::String, filename::String, axis_labels::Vector{String}, data::Vector{V}; kwargs...) where V <: AbstractVector

Saves data to a CSV file with metadata and headers. The first line of the CSV contains serialized metadata, the second line contains column headers, and subsequent lines contain the data.

# Arguments
- `typeofdata::Symbol`: Specifies the type of data (e.g., `:ℝ_to_ℝ`, `:ℝ³_to_ℝ`, or `"bandstructure"`).
- `path::String`: Directory path where the file will be saved.
- `filename::String`: Name of the CSV file.
- `axis_labels::Vector{String}`: Labels for the data columns (e.g., axis names).
- `data::Vector{V}`: Data to be saved, where each element is a vector corresponding to a column.
- `kwargs...`: Optional keyword arguments. Supported:
    - `row_input`: If `true`, preprocesses data as row input.

# Metadata Format
The metadata (first line of the CSV) is a serialized and base64-encoded Julia object containing:
- Type of plot/data
- Number of entries
- Codomain dimensionality
- Additional keyword arguments

# Data Format
- Second line: Column headers (from `axis_labels`)
- Third line onwards: Data rows

# Supported `typeofdata` Values
- `:ℝ_to_ℝ`: 1D to 1D data (e.g., x and y columns)
- `:ℝ³_to_ℝ`: 3D to 1D data (e.g., x, y, z, and C columns)
- `"bandstructure"`: Special handling for band structure data

# Example
"""
function save_data(typeofdata::Symbol, path::String, filename::String, axis_labels::Vector{String}, data::Vector{V}; kwargs...) where V <: AbstractVector
    # this code will save data in the following format, where the first line of a CSV is the metadata
    # the second line would be the headers for the data
    # and the third row onwards would be data
    # Ex:
    # |typeofdata::String (Ex: ℝ_to_ℝ)| npts::Int | dims_codomain::Int (n as in ℝ→ℝⁿ) | title_of_plot::String | flip x and y axes?
    # | xlabel  | ylabel|zlabel | etc | etc|
    # |data     | data  | data  | data|data|

    row_input = get(kwargs, :row_input, false)
    
    if row_input
        preproc = reduce(vcat,transpose.(positions))
        data = [preproc[:,i] for i = 1:size(preproc)[2]]
    end
    if typeofdata == :ℝ_to_ℝ
        # metadata is in the format ["type of plot", n_entries, n_codomain (dimensionality of codomain)]
        metadata = [Symbol(typeofdata), size(data[1])[1], 1, kwargs]
        df = DataFrame(x = data[1], y = data[2])
    elseif typeofdata == :ℝ³_to_ℝ
        x = data[1]
        y = data[2]
        z = data[3]
        C = data[4]
        # metadata is in the format ["type of plot", n_entries, n_codomain (dimensionality of codomain)]
        metadata = [Symbol(typeofdata), size(x)[1], 1, kwargs]
        df = DataFrame(x = x, y = y, z=z, C=C)
    elseif typeofdata == "bandstructure"
        # total number of kpts that hamiltonian has been sampled over
        nkpts = size(data[2])[1]
        # total number of labels that will go on the x axis
        nklabels = size(klabels)[1] 
        # number of interpolation points between kpts
        ninterpolate = Int((nkpts-1)/(nklabels-1))
        kxs = Float64[]
        kys = Float64[]
        kzs = Float64[]
        # okay, we will now make column 1 which will hold the k point list
        bands_labels = String[] 
        for ik = 1:(nklabels-1)
            push!(bands_labels,klabel[ik])
            for ik = 2:nk
                push!(bands_labels,"")
            end
        end
        push!(bands_labels,klabel[end])
        # make columns 2-4
        for k ∈ data[2] #loop over all of the k vectors
            push!(kxs,k[1])
            push!(kys,k[2])
            push!(kzs,k[3])
        end
        Evals = data[3]
        nbands = size(Evals)[2]
        # now make the labels and columns 5-N
        labels = ["k point", "kx (1/m)", "ky (1/m)", "kz (1/m)"]
        for iE = 1:nbands
            push!(labels,"E_"*string(iE))
        end
        # okay, now figure out the code to get this stuff into the data DataFrame
        # each column for the energy_values = Evals[ik, :]
        # each row should correspond to one k-point 
        
    end
    if(@isdefined df)
        rename!(df, Symbol.(axis_labels))
        save_csv(path, filename, df, metadata)
    else
        println("Saving failed")
    end
end

function save_csv(output_dir::String, filename::String, df::DataFrame, metadata::Vector{Any})
    # Serialize to an IOBuffer and write base64 string to file
    buffer = IOBuffer()
    serialize(buffer, metadata)
    serialized_data = take!(buffer)
    base64_string = base64encode(serialized_data)

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
            println(file, base64_string)
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
