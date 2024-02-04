using CSV
using DataFrames


function get_data(readDir::String)
    println("Reading from: ", readDir)
    try
        csv_file_path = joinpath(readDir, "scatterplot.csv")
        df = CSV.File(csv_file_path) |> DataFrame
        # println(df)
        x = df.X
        y = df.Y
        z = df.Z
        d = df.D

        # plt = scatter(x, y, z, color = d, marker_z = d, legend = false)
        # display(plt)
        # println(d)
        return x, y, z, d
    catch e
        println("An IO error probably occurred: ", e)
        return nothing
    end
end

# get_data()