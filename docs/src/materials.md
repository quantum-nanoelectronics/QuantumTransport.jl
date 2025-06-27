# materials.jl

## Materials

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

### Metal Hopping

This function calculates the hopping terms for metallic materials within a lattice structure and adds them to a vector of Hopping objects.

#### Tutorial

To compute hopping terms for metallic materials within a lattice structure, follow these steps:

1. Iterate over neighboring lattice sites to compute hopping energies.
2. Compute the orbital index `iorb` from the initial lattice site indices `ia`.
3. Compute the hopping energy `t`.
4. Add the hopping term corresponding to the on-site energy between the initial lattice site `ia` and itself to the vector `NNs`.
5. Iterate over all three spatial dimensions and both positive and negative directions.
6. Compute the shift vector `di` to determine neighboring lattice sites.
7. Compute the hopping energy `t` and add the hopping term between the initial lattice site `ia` and the neighboring lattice site `ib` to the vector `NNs`.

**Function**

```julia
metalHopping(p::NamedTuple, NNs::Vector{Hopping}, ia::Vector{Int})
```

---

### insHopping

This function calculates the hopping terms for insulating materials within a lattice structure and adds them to a vector of Hopping objects.

#### Tutorial

1. Define the parameters for the material using a named tuple `p`.
2. Create an empty vector `NNs` to store the hopping terms.

```julia
NNs = Vector{Hopping}()
```

3. Call the `insHopping` function with the parameters `p`, `NNs`, and the initial lattice site indices `ia`.

```julia
insHopping(p, NNs, ia)
```

4. Iterate over neighboring lattice sites to compute hopping energies.
5. Compute the orbital index `iorb` from the initial lattice site indices `ia`.
6. Compute the hopping energy `t` using the parameters `p.ϵ₁` and `p.t`.
7. Add the hopping term corresponding to the on-site energy between the initial lattice site `ia` and itself to the vector `NNs`.
8. Iterate over all three spatial dimensions and both positive and negative directions.
9. Compute the shift vector `di` to determine neighboring lattice sites.
10. Compute the hopping energy `t` and add the hopping term between the initial lattice site `ia` and the neighboring lattice site `ib` to the vector `NNs`.

---

### pushHopping

Appends a hopping term to a vector of hopping terms (`NNs`).

#### Tutorial

1. Convert the initial and final lattice site indices (`ia` and `ib`) to their corresponding orbital indices (`a` and `b`) using the function `xyztoi`.
2. Convert the initial and final lattice site indices (`ia` and `ib`) to their corresponding real-space coordinates (`ra` and `rb`) using the function `xyztor`.
3. Calculate the displacement vector `r` between the initial and final sites by subtracting their real-space coordinates (`rb - ra`).
4. Create a new `Hopping` object containing the hopping information.
5. Append the newly created `Hopping` object to the vector `NNs` using the `push!` function.

---

### weyl3Hopping

Calculates hopping terms for a 3D material with Weyl points in its band structure, considering nearest, next-nearest, and next-to-next-nearest neighbors.

#### Tutorial

1. Define initial orbital and lattice site indices.
2. **Nearest neighbors**: Iterate and compute hopping tensors using Pauli matrices based on axis and direction.
3. **Next-nearest neighbors**: Use combinations of displacements in 2D to define hopping.
4. **Next-to-next-nearest neighbors**: Use doubled displacements and compute hopping energy using Pauli matrices.

---

### weyl2Hopping

Calculates hopping terms for a 2D material with Weyl points in its band structure, considering nearest and next-nearest neighbors.

#### Tutorial

1. Define initial orbital and lattice site indices.
2. **Nearest neighbors**: Iterate and compute hopping tensors based on axis and direction.
3. **Next-nearest neighbors**: Use combinations of displacements to define hopping energy using Pauli matrices.

---

### chern2DHopping

Calculates hopping terms for a 2D material with Chern points, considering on-site energy, nearest-neighbor coupling mediated by spin-orbit interaction, and normal hopping between nearest neighbor orbitals.

#### Tutorial

1. **On-site term**: Use `nextsite(iorb)` and material parameters to compute `t`.
2. **Nearest-neighbor term**: Use displacement vectors and compute `t` with spin-orbit interaction using Pauli matrices.
3. **Normal hopping term**: Restore original orbital index and compute direct hopping energy using identity matrix.