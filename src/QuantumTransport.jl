module QuantumTransport

# Export the function if you want it to be accessible from outside the module, such as the test folder
# Include the file where the function is defined

# Run sample test codes
include("sample-code/hello-world.jl")
include("sample-code/column-major.jl")
export sampleTest, calculate_time, row_major, column_major


# Run matrix inversion test codes
include("matrices/woodbury-inverse-1.jl")
include("matrices/recursive-greens/RGF.jl")
export block_inv_main, rgf_main

# Run input output test codes
include("io/read-positions.jl")
include("io/write-positions.jl")
export generate_csv, get_data

# Run data visualization test codes
include("data-vis/plot.jl")
export generate_plot_makie, plot_pos

end
