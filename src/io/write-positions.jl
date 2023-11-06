using CSV
using DataFrames
using Random

include("../data-vis/make-cnt-data.jl")

const IO_DIR = "src/io/"

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



function generate_csv(output_dir::String=@__DIR__, positions::Vector{Vector{Float64}}=SAMPLE_POSITIONS)
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

# generate_csv(IO_DIR, positions)

# generate_csv(IO_DIR)




###
### Archive
###

# Random.seed!(123)
# num_rows = 100
# num_cols = 3
# dataRand = Vector{Vector{Float64}}(undef, num_rows)

# for i in 1:num_rows
#     row = [rand(1.0:100.0) for _ in 1:num_cols]
#     push!(dataRand, row)
# end


# generate_csv(dataRand)
# function generate_csv(positions::Vector{Vector{Float64}})
#     csv_file_path = "positions.csv"
#     data = []
#     df = DataFrame(positions, [:X, :Y, :Z])
#     println(df)

#     for i = 1:length(positions[1])
#         # x, y, z = positions[i]
#         # randomElectronDensity = x + y + z
#         # println("Position $x $y $z | Density $randomElectronDensity")
#         randomElectronDensity = rand(1.0:100.0)
#         push!(data, randomElectronDensity)
        
#     end

#     df[:, :D] = data
#     # df = DataFrame(Float64(data))
#     #println(data[:, 1])
#     #column_names = [:x, :y, :z, :randomElectronDensity]  # Specify column names
#     #df = DataFrame()  # Create a DataFrame with mixed data types

#     csv_file_path = "scatterplot.csv"
#     CSV.write(csv_file_path, df)
#     println("CSV file '$csv_file_path' has been created.")
# end



