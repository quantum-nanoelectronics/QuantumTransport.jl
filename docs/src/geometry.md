# geometry.jl

This module provides functions for determining the material type based on position coordinates.

## Tutorial

### Example: Material Determination

To determine the material type based on position coordinates, follow these steps:

1. **Define the position coordinates using a `Vector{Float64}`:**
    ```julia
    R = [x, y, z]
    ```

2. **Call the `geometry` function to determine the material type:**
    ```julia
    material_type = geometry(R)
    ```
    The `geometry` function returns a string indicating the material type (`"insulator"` or `"GaAs"`).

---

## Constants

- `geometry_params`: Named tuple containing parameters for defining the geometry of the device.
    - `A`: Matrix of unit cell lattice vectors.
    - `nx`: Number of times to tile the cell over space in the x-direction.
    - `ny`: Number of times to tile the cell over space in the y-direction.
    - `nz`: Number of times to tile the cell over space in the z-direction.
    - `prune`: List of dimensions to prune.

### Example Implementation

```julia
function geometry(R::Vector{Float64})
    x = R[1]; y = R[2]; z = R[3];
    if (x < 10*nm || x > 50*nm)
        return "insulator"
    else
        return "GaAs"
    end
end

geometry_params = (
    A = 2.866 * nm * I(3),   # Matrix of unit cell lattice vectors
    nx = 50,                 # Number of times to tile the cell over space in the x-direction
    ny = 5,                  # Number of times to tile the cell over space in the y-direction
    nz = 1,                  # Number of times to tile the cell over space in the z-direction
    prune = ["x", "y", "z"], # List of dimensions to prune
)
```