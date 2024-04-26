# Data visualization module
# Adding Makie, etc. to the dependencies will not work on github
#using Makie #If this line is uncommented or if Makie is added to this package, github tests will fail
module DataVisualizationModule

using ..CommonModule
using ..InputOutputModule

using CSV
using DataFrames
# using Plots
# using Makie # Cannot do this 
#using Makie #If this line is uncommented or if Makie is added to this package, github tests will fail

using ColorSchemes
using Colors

using GLMakie
using CairoMakie
# Adding Makie, etc. to the dependencies will not work on github
#using Makie #If this line is uncommented or if Makie is added to this package, github tests will fail

#using Makie #If this line is uncommented or if Makie is added to this package, github tests will fail
include("plot.jl")
include("plotstuff.jl")
export generate_plot_makie, plot_pos, call_function_based_on_header

end