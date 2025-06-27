# runs.jl

Defines parameters for different types of simulations.

### Tutorial: Running Simulations

1. Set the output path:
   ```julia
   path = "./"
   ```

2. Define simulation modes:

#### ➤ Unit Cell Simulation
```julia
unitcell = (
    bands = true,
    bands_project = [σ[1], γ⁵],
    poisson = false,
    DOS = false
)
```

#### ➤ Transport Simulation
```julia
transport = (
    ΔV = 0.05,
    μ = 0.1 * eV,
    T = 300,
    η = 1e-4 * eV,
    savedensities = true,
    density_project = [I(2), [σ[1], σ[2], σ[3]]],
    Gʳinv_method = "CBR",
    D_dephasing = 0.1 * eV,
    D_spin = 0.01 * eV,
    D_momentum = 0.5 * eV,
    kspace = false
)
```

#### ➤ Supercell Simulation
```julia
supercell = (
    bands = true,
    bands_project = [σ[1], σ[2]],
    poisson = true,
    μ = 0.1 * eV,
    T = 300,
    η = 1e-4 * eV,
    savedensities = true,
    density_project = [I(2), [σ[1], σ[2], σ[3]]],
    Gʳinv_method = "CBR",
    D_dephasing = 0.1 * eV,
    D_spin = 0.01 * eV,
    D_momentum = 0.5 * eV
)
```

3. Create a named tuple with the simulation setup:
```julia
runparams = (
    path = path,
    unitcell = unitcell,
    transport = transport,
    supercell = supercell
)
```

4. Run the simulation:
```julia
run_simulation(runparams)
```