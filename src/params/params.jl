# set ϵ₁ in params with ε₁ 
const hoppings = ["n","nx","ny","nz","nsite","norb","t","SLa₂","a₁","a₂","a₃",
"deviceMaterial","ε₁","A","fieldtype"]
const hamiltonian = ["μ","μ_disorder", "n", "nx", "ny", "nz", "norb", "nsite", 
"A", "t", "ε₁", "a₁","a₂","a₃", "deviceMagnetization", "fieldtype"]
const hoppingMatrix = ["n", "nsite", "norb", "t", "ε₁"]
const electrodes = ["electrodeMaterial","nsite","norb","η"]

# ϵ₁
# t₉ = 0.0, 
# vf = 1000000, 
# η = 0.0001, 
# ε = 0.0, 
# ϵ₁ = 2.0, 
function generateParams(params::Dict)
    nested_params = Dict(Dict())
    hopping_params = Dict()
    hamiltonian_params = Dict()
    matrix_params = Dict()
    electrode_params = Dict()

    for key in hoppings
        hopping_params[key] = params[key]
    end

    for key in hamiltonian
        hamiltonian_params[key] = params[key]
    end

    for key in hoppingMatrix
        matrix_params[key] = params[key]
    end

    for key in electrodes
        electrode_params[key] = params[key]
    end

    nested_params["hoppings"] = hopping_params
    nested_params["hamiltonian"] = hamiltonian_params
    nested_params["matrix"] = matrix_params
    nested_params["electrode"] = electrode_params
    return nested_params
end