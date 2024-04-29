module TestDataVisualization
using QuantumTransport
using Test
#using Makie #If this line is uncommented or if Makie is added to this package, github tests will fail

function runVisualizationTests()
    readDir = OUTPUT_DIR 
    filename = "transmission.csv"  # specify the CSV file name

    file_path = joinpath(readDir, filename)
    fig = call_function_based_on_header(readDir, filename)
    # filename = "CNTpositions.csv"
    # call_function_based_on_header(readDir, filename)
    # Save the figure
    @test isfile(file_path)
    return fig
end

function runVisualizationTestsGLMakie()
    readDir = OUTPUT_DIR 
    filename = "transmission.csv"  # specify the CSV file name

    file_path = joinpath(readDir, filename)
    call_function_based_on_header(readDir, filename) # pass in true to display with GLMakie backend
    filename = "CNTpositions.csv"
    call_function_based_on_header(readDir, filename)
    # Save the figure
    @test isfile(file_path)
end

runVisualizationTests()

end # module TestDataVisualization
