# Reading csv files:
# The first row of the csv file should contain information about the type of plot.
# The second row of the file should contain the column headings.
# The third row onwards should contain the data.

module InputOutputModule # the _ in the directory is needed because this needs to be compiled first to be used within other package modules

include("readfiles.jl")
include("writefiles.jl")
export get_data, save_csv, save_data

end