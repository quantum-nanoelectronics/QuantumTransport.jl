# Data visualization module
# Adding Makie, etc. to the dependencies will not work on github
#using Makie #If this line is uncommented or if Makie is added to this package, github tests will fail
module DataVisualizationModule

using ..QuantumTransport
using CSV
using DataFrames
using ColorSchemes
using Colors
using CairoMakie
using Dates
# Adding Makie, etc. to the dependencies will not work on github
#using Makie #If this line is uncommented or if Makie is added to this package, github tests will fail

include("plot.jl")
include("functions.jl")
export plot, get_empty_plot, save_plot

end