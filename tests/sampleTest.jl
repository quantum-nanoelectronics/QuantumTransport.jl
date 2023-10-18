using Test
include("../src/plotting/testPlot.jl")

@testset "Test generate_plot_makie function" begin
    @test_throws Exception generate_plot_makie(-1)  # Testing with a negative margin
    @test_throws Exception generate_plot_makie(0)   # Testing with zero margin
    @test generate_plot_makie(0.1) !== nothing      # Testing with a valid margin
end
