# Data visualization module

module DataVisualizationModule

using ..CommonModule

#using Makie #If this line is uncommented or if Makie is added to this package, github tests will fail
include("plot.jl")
include("plotstuff.jl")
export generate_plot_makie, plot_pos, call_function_based_on_header

end