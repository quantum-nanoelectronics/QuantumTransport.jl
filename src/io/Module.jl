# reading csv files:
# The first row of the file should contain the column headings.
# The second row of the csv file should contain information about the type of plot.
# The third row onwards should contain the data.

module InputOutputModule

include("read-positions.jl")
include("write-positions.jl")
export generate_csv, get_data, save_csv

end