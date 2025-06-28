module TestScattering

using QuantumTransport
using Test

# Bypassing Driver for a sweep accross ϵ_rand_strengths

include(joinpath(INPUT_DIR, "AllInputs.jl"))

p = runparams["transport"]
ϵ_values = [0.0, 0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]

for ϵ in ϵ_values
    p["ϵ_rand_strength"] = ϵ

    # Run transport
    transport(p)

    # Visualize this transport run
    include("TestVisualization.jl")

end

end # module