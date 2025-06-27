# Geometry Documentation

This documentation describes the input file used to define device geometry, simulation parameters, and material properties in a quantum transport simulation framework.

---

## Module: `geometry.jl`

Provides functions to determine material type based on position coordinates.

### Tutorial: Material Determination

To determine the material type at a given point:

1. Define position coordinates:
   ```julia
   R = [x, y, z]  # Vector{Float64}
   ```

2. Call the geometry function:
   ```julia
   material_type = geometry(R)
   ```

3. The result is a `String`: either `"insulator"` or `"GaAs"`.

### Constants

- `geometry_params`: Named tuple of device geometry settings:
  - `A`: Unit cell lattice vector matrix
  - `nx`, `ny`, `nz`: Number of tiled unit cells in each dimension
  - `prune`: Dimensions to prune

### Example Code

```julia
module geometry

using FANCY_TRANSPORT_PACKAGE
export geometry, geometry_params

function geometry(R::Vector{Float64})
    x, y, z = R
    if (x < 10 * nm || x > 50 * nm)
        return "insulator"
    else
        return "GaAs"
    end
end

geometry_params = (
    A = 2.866 * nm * I(3),
    nx = 50,
    ny = 5,
    nz = 1,
    prune = ["x", "y", "z"],
)

end  # module
```

---

## Module: `runs.jl`

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

---

## Module: `materials.jl`

Provides functionality to compute material-specific hopping terms.

### Example Usage

1. Define material parameters:
   ```julia
   p = ( ... )  # Named tuple
   ```

2. Create a container for hopping terms:
   ```julia
   NNs = Vector{Hopping}()
   ```

3. Call the material-specific hopping function:
   ```julia
   metalHopping(p, NNs, ia)
   ```

- `p`: Named tuple of material parameters
- `NNs`: Vector of `Hopping` objects
- `ia`: Indices of the initial lattice site

---