# QuantumTransport.jl

QuantumTransport.jl is a Julia package designed for simulating quantum transport phenomena in nanoscale electronic devices.

## Status

[![CI Status](https://github.com/quantum-nanoelectronics/QuantumTransport.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/quantum-nanoelectronics/QuantumTransport.jl/actions)


This package supports Julia 1.9.3 and above on Linux, MacOS, and Windows.

## Installation

To install QuantumTransport.jl, open the Julia REPL, type `]` to enter package mode, and run:

```julia
pkg> add QuantumTransport
```

Alternatively, using the `Pkg` API:

```julia
julia> import Pkg

julia> Pkg.add("QuantumTransport")
```

## Usage

After [installing](#installation) QuantumTransport, you can start using the package with:

```julia
using QuantumTransport

# Your simulation code goes here
```

See the `examples/` directory for real-world simulations and usage patterns.

## Contributing and Support

We welcome contributions! Please open a new [pull request](https://github.com/quantum-nanoelectronics/QuantumTransport.jl/pulls) from a feature branch. 

If you have a bug report, feature request, or general question, please open an [issue](https://github.com/quantum-nanoelectronics/QuantumTransport.jl/issues).

### Running Tests

If you're developing the package locally, you can run the test suite from the Julia Pkg REPL:


```julia
pkg> activate /your/path/to/the/cloned/QuantumTransport.jl/.

pkg> instantiate

pkg> test
```
