module TestDataVisualization
using QuantumTransport
using Test
using CairoMakie

"""
    runVisualizationTests()

Run the visualization tests for the QuantumTransport package.
Cannot fully test Data Visualization on GitHub, as no GPU is available.
Non-interactive image plots generated, for now.

"""
function runVisualizationTests()
    nm = 1E-9
    margin = 0.0
    filename = "scatterplot.csv"
    baseDir = abspath(joinpath(@__DIR__, ".."))
    ioDir = joinpath(baseDir, "data-output")
    vals = get_data(ioDir, filename)
    df = vals[1]
    metaData = vals[2]
    
    println("-Reading-")
    println("DataFrame: ")
    println(first(df, 5))
    println("Metadata: ")
    println(metaData)
    
    
    # Obtain the figure object
    figure = generate_plot_makie(df, margin)
    @test figure isa Figure
    # Save the figure
    file_path = joinpath(ioDir, "figure_from_generate_plot_makie.png")
    save(file_path, figure)
    @test isfile(file_path)
    
    # Obtain the figure object
    figure = plot_pos(df, nm)
    @test !isnothing(figure)
    # Save the figure
    file_path = joinpath(ioDir, "figure_from_plot_pos.png")
    save(file_path, figure)
    @test isfile(file_path)
    
end

runVisualizationTests()

end # module TestDataVisualization
