using CSV
using DataFrames

# Returns the DataFrame and metadata from the scatterplot.csv file in the given directory
function get_data(readDir::String)
    println("Reading from: ", readDir)
    try
        csv_file_path = joinpath(readDir, "scatterplot.csv")
        df = CSV.File(csv_file_path) |> DataFrame

        metadata = collect(df[1, :])
        delete!(df, 1)

        metadata_array = [String(x) for x in metadata if x !== missing]

        return df, metadata_array
    catch e
        println("An IO error probably occurred: ", e)
        return nothing, nothing
    end
end

# get_data()