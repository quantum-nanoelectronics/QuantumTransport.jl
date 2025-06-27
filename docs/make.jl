using Documenter
using QuantumTransport

Documenter.DocMeta.setdocmeta!(QuantumTransport, :DocTestSetup, :(using QuantumTransport); recursive=true)

makedocs(
    sitename = "QuantumTransport.jl",
    modules = [QuantumTransport],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "geometry.jl" => "geometry.md",
        "materials.jl" => "materials.md",
        "runs.jl" => "runs.md",
    ],
)

deploydocs(
    repo = "github.com/quantum-nanoelectronics/QuantumTransport.jl",
    push_preview = true,
)