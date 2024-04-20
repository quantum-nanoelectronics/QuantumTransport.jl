module TestDataVisualization
using QuantumTransport
using Test

function runVisualizationTests()
    
    readDir = OUTPUT_DIR 
    filename = "transmission.csv"  # specify the CSV file name

    file_path = joinpath(readDir, filename)
    call_function_based_on_header(readDir, filename)
    filename = "CNTpositions.csv"
    call_function_based_on_header(readDir, filename)
    # Save the figure
    @test isfile(file_path)
    
    
end

runVisualizationTests()

end # module TestDataVisualization
