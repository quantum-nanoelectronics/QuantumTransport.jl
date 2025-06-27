# geometry.jl

This documentation describes the input file used to define device geometry, simulation parameters, and material properties in a quantum transport simulation framework.

---

## Input File: `geometry.jl`

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

using QuantumTransport
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
