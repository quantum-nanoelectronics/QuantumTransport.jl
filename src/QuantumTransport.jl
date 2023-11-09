module QuantumTransport

# Export the function if you want it to be accessible from outside the module, such as the test folder
# Include the file where the function is defined

# Hello World function
export sampleTest
include("hello-world.jl")

# Column major arrays
export calculate_time, row_major, column_major
include("column-major.jl")

export block_inv_main
include("matrices/woodbury-inverse-1.jl")

export generate_csv, get_data
include("io/read-positions.jl")
include("io/write-positions.jl")

export generate_plot_makie, plot_pos
include("data-vis/plot.jl")

end
