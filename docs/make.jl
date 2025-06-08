using Documenter
using QuantumTransport

makedocs(
    sitename = "QuantumTransport.jl",
    modules = [QuantumTransport],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "geometry.jl" => "geometry.md",
        "materials.jl" => "materials.md",
    ],
)

deploydocs(
    repo = "github.com/quantum-nanoelectronics/QuantumTransport.jl",
    push_preview = true,
)