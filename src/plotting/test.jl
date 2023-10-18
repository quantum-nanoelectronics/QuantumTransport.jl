using Pkg
#Pkg.add("CSV")
#Pkg.add("DataFrames")
#Pkg.add("Plots")

#using Plots
using CSV
using DataFrames
using Random




function generate_csv(positions::Vector{Vector{Float64}})
    csv_file_path = "positions.csv"
    data = []
    df = DataFrame(positions, [:X, :Y, :Z])
    println(df)

    for i = 1:length(positions[1])
        # x, y, z = positions[i]
        # randomElectronDensity = x + y + z
        # println("Position $x $y $z | Density $randomElectronDensity")
        randomElectronDensity = rand(1.0:100.0)
        push!(data, randomElectronDensity)
        
    end

    df[:, :D] = data
    # df = DataFrame(Float64(data))
    #println(data[:, 1])
    #column_names = [:x, :y, :z, :randomElectronDensity]  # Specify column names
    #df = DataFrame()  # Create a DataFrame with mixed data types

    csv_file_path = "scatterplot.csv"
    CSV.write(csv_file_path, df)
    println("CSV file '$csv_file_path' has been created.")
end



data = [
    [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], # X
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], # Y
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], # Z
    ]
generate_csv(data)


# Random.seed!(123)
# num_rows = 100
# num_cols = 3
# dataRand = Vector{Vector{Float64}}(undef, num_rows)

# for i in 1:num_rows
#     row = [rand(1.0:100.0) for _ in 1:num_cols]
#     push!(dataRand, row)
# end


# generate_csv(dataRand)




