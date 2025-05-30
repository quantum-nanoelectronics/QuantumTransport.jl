# QuantumTransport.jl

QuantumTransport.jl is a Julia package designed for simulating quantum transport in nano-electronic devices.

## Status

CI Status: 

[![Build Status](https://github.com/quantum-nanoelectronics/QuantumTransport.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/quantum-nanoelectronics/QuantumTransport.jl/actions/workflows/CI.yml?query=branch%3Amain)

View [more](https://github.com/quantum-nanoelectronics/QuantumTransport.jl/actions) CI on GitHub Actions


## Installation

To install QuantumTransport from the Julia's REPL, type `]` to enter the Pkg REPL and run:

```julia
pkg> add QuantumTransport
```

Or, equivalently:

```julia
julia> import Pkg

julia> Pkg.add("QuantumTransport")
```

## Usage

After installing QuantumTransport, you can start using the package by:

```julia
using QuantumTransport

# Your simulation code goes here
```

For detailed examples and usage instructions, please refer to the `examples` directory.

## Questions and Contributions

Contributions are welcome, as are feature requests! If you'd like to contribute, please open a pull request and use a feature branch. Please open an [issue](https://github.com/quantum-nanoelectronics/QuantumTransport.jl/issues) if you encounter any problems with the package. 

### Running Tests for Developers

To test QuantumTransport from the Julia's REPL after cloning this repository, type `]` to enter the Pkg REPL and run:

```julia
pkg> activate /path/to/your/QuantumTransport.jl/

pkg> instantiate

pkg> test
```

This package was built on Julia 1.9.3. 
