using CSV
using DataFrames
using Base.Filesystem: mktemp

function get_data(readDir::String)
    println("Reading from: ", readDir)
    csv_file_path = joinpath(readDir, "scatterplot.csv")
    
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

