using Documenter
using QuantumTransport

makedocs(
    sitename = "QuantumTransport.jl",
    modules = [QuantumTransport],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
    ],
)

deploydocs(
    repo = "github.com/quantum-nanoelectronics/QuantumTransport.jl",
    push_preview = true,
)