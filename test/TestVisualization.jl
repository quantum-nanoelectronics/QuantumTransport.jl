module TestDataVisualization
using QuantumTransport
using Test
#using Makie #If this line is uncommented or if Makie is added to this package, github tests will fail

function runVisualizationTests()
    readDir = OUTPUT_DIR 
    filename = first(filter(f -> startswith(f, "transmission"), readdir(readDir)))
    println("filename= ", filename)
    file_path = joinpath(readDir, filename)
    fig = plot(readDir, filename)
    # filename = "CNTpositions.csv"
    # plot(readDir, filename)
    # Save the figure
    @test isfile(file_path)
    return fig
end

function runVisualizationTestsGLMakie()
    readDir = OUTPUT_DIR 
    filename = "transmission.csv"  # specify the CSV file name

    file_path = joinpath(readDir, filename)
    plot(readDir, filename) # pass in true to display with GLMakie backend
    filename = "CNTpositions.csv"
    plot(readDir, filename)
    # Save the figure
    @test isfile(file_path)
end


function recordVisualization()
    readDir = OUTPUT_DIR
    filename = "CNTpositions.csv"
    fig, ax, plt = plot(readDir, filename)
    
    record(fig, "visualization.mp4", 1:120) do frame
        ax.azimuth[] = 1.7pi + 0.3 * sin(2pi * frame / 120) # Camera tilt effect
        # Update other necessary components if required
    end
end

runVisualizationTests()

end # module TestDataVisualization
