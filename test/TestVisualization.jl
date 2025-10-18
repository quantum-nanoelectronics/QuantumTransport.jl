module TestDataVisualization
using QuantumTransport
using Test
#using Makie #If this line is uncommented or if Makie is added to this package, github tests will fail

include("TestVisualizationFunctions.jl")

runVisualizationTests()

end # module TestDataVisualization
