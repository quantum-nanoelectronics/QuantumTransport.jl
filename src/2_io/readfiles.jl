using CSV
using DataFrames
using Base.Filesystem: mktemp

"""
    get_data(readDir::String, filename::String)

Reads data from a file.
The first row of the csv file should contain information about the type of plot.
The second row of the file should contain the column headings.
The third row onwards should contain the data.

# Arguments
- `readDir::String`: The directory path where the file is located.
- `filename::String`: The name of the file to read.

# Returns
The data read from the file.

"""
function get_data(readDir::String, filename::String)
    csv_file_path = joinpath(readDir, filename)
    
    # Open the original CSV file to read the metadata line
    metadata_line = open(csv_file_path, "r") do file
        readline(file)  # This reads the first line, which is the metadata
    end
    # Split the metadata into an array
    metadata_array = split(metadata_line, ",")

    # Create a temporary file to write the CSV content excluding the first metadata line
    temp_csv_path, temp_file = mktemp()
    close(temp_file)  # Close the temp file handle since it's automatically opened by mktemp
    
    # Read the original file and write to the temporary file, skipping the first line
    open(csv_file_path, "r") do infile
        readline(infile)  # Skip the metadata line
        content = read(infile)  # Read the rest of the file
        open(temp_csv_path, "w") do outfile
            write(outfile, content)  # Write the content to the temp file
        end
    end

    # Now read the DataFrame from the temporary file
    df = CSV.File(temp_csv_path) |> DataFrame

    # Clean up: remove the temporary file
    rm(temp_csv_path; force=true)

    return df, metadata_array
end

