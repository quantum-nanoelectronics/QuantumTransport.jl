# runs.jl

This file contains parameters for configuring different types of simulations.

---

## Tutorial

### Example: Running Simulations

To run simulations with different configurations, follow these steps:

1. **Define the directory path where simulation data will be stored:**
    ```julia
    path = "./"
    ```

2. **Configure parameters for simulations focusing on the electronic properties of a single unit cell:**
    ```julia
    unitcell = (
        bands = true,
        bands_project = [σ[1], γ⁵],
        poisson = false,
        DOS = false
    )
    ```
    - `bands`: Boolean indicating whether to compute electronic band structure.
    - `bands_project`: List of projection operators for computing bands.
    - `poisson`: Boolean indicating whether to include Poisson solver.
    - `DOS`: Boolean indicating whether to compute density of states.

3. **Configure parameters for simulations focusing on voltage-dependent transport properties:**
    ```julia
    transport = (
        ΔV = 0.05,
        μ = 0.1 * eV,
        T = 300,
        η = 1E-4 * eV,
        savedensities = true,
        density_project = [I(2), [σ[1], σ[2], σ[3]]],
        Gʳinv_method = "CBR",
        D_dephasing = 0.1 * eV,
        D_spin = 0.01 * eV,
        D_momentum = 0.5 * eV,
        kspace = false
    )
    ```
    - `ΔV`: Voltage bias.
    - `μ`: Chemical potential.
    - `T`: Temperature.
    - `η`: Broadening parameter for densities.
    - `savedensities`: Boolean indicating whether to save density matrices.
    - `density_project`: List of projection operators for density matrices.
    - `Gʳinv_method`: Method for computing Green's function inverse.
    - `D_dephasing`: Dephasing parameter.
    - `D_spin`: Spin relaxation parameter.
    - `D_momentum`: Momentum relaxation parameter.
    - `kspace`: Boolean indicating whether to use k-space.

4. **Configure parameters for simulations using multiple unit cells:**
    ```julia
    supercell = (
        bands = true,
        bands_project = [σ[1], σ[2]],
        poisson = true,
        μ = 0.1 * eV,
        T = 300,
        η = 1E-4 * eV,
        savedensities = true,
        density_project = [I(2), [σ[1], σ[2], σ[3]]],
        Gʳinv_method = "CBR",
        D_dephasing = 0.1 * eV,
        D_spin = 0.01 * eV,
        D_momentum = 0.5 * eV
    )
    ```
    (Same parameters as unitcell and transport)

5. **Create a named tuple `runparams` containing the configured parameters:**
    ```julia
    runparams = (path = path, unitcell = unitcell, transport = transport, supercell = supercell)
    ```

6. **Use the `run_simulation` function with the `runparams` tuple to execute the simulation:**
    ```julia
    run_simulation(runparams)
    ```

---

## Constants

- `path`: Path to the directory where simulation data will be stored.
