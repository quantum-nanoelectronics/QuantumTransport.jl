module ElectrodesModule

export Electrode, genΣₖs

include("../hoppings/createHoppings.jl")
include("Utilities.jl")
include("ModifyH.jl")
include("SelfEnergies.jl")


# from Constants.jl
ħ = 1.05457E-34
h = ħ * 2*π
m₀ = 9.10938E-31
q = 1.60218E-19
ϵ₀ = 8.854E-12
metre = 1
au = 27.2
eV = 1.0
cm = 1E-2
μ₀ = 1.2566*10^-6 # H/m vacuum permeability
nm = 1E-9
Å = 1E-10
r₀ = 5.29E-11
Ry = m₀*q^4 / (8*h^2*ϵ₀^2)
μₑ = 9.28*10^-24 # electron magnetic moment in A*m^2
μB = 5.788838E-5 # bohr magneton in eV/T
kB = 8.617E-5 # boltzmann constant in eV/K

end