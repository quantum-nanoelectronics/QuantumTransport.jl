# reading csv files:
# The first row of the csv file should contain information about the type of plot.
# The second row of the file should contain the column headings.
# The third row onwards should contain the data.

module InputOutputModule

include("read-positions.jl")
include("write-positions.jl")
export get_data, save_csv

end