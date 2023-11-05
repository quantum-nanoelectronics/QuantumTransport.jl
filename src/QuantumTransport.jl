module QuantumTransport

# Export the function if you want it to be accessible from outside the module, such as the test folder

# Hello World function
export sampleTest
# Include the file where the function is defined
include("hello-world.jl")

# Column major arrays
export calculate_time, row_major, column_major
include("column-major.jl")


end
