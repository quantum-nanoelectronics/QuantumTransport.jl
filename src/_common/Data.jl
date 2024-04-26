BASE_DIR = abspath(joinpath(@__DIR__, "../.."))
INPUT_DIR = joinpath(BASE_DIR, "data-input")
OUTPUT_DIR = joinpath(BASE_DIR, "data-output")

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

#unit conversions
μₑ = 9.28*10^-24 # electron magnetic moment in A*m^2
#Å = 1.8897 * r0 	#angstrom, in hartree units
nm = 10 * Å 		#angstrom, in hartree units
kT_RT = 0.02585*eV
#metre = 10^10 * Å


τ₀ = [
    1 0
    0 1
]

τ₁ = [
    0 1
    1 0
]

τ₂ = [
    0 -im
    im 0
]

τ₃ = [
    1 0
    0 -1
]

σ₀ = [
    1 0
    0 1
]

σ₁ = [
    0 1
    1 0
]
σ₂ = [
    0 -im
    im 0
]
σ₃ = [
    1 0
    0 -1
]


S1z = [1 0 0;
    0 0 0;
    0 0 -1]

S1y = sqrt(1 / 2) * [0 -im 0;
    im 0 -im;
    0 im 0]

S1x = sqrt(1 / 2) * [0 1 0;
    1 0 1;
    0 1 0]


S² = [
    (1/2)*(1/2+1) 0
    0 (1/2)*(1/2+1)
]


σ = Vector{Matrix}(undef, 3)
σ[1] = σ₁;
σ[2] = σ₂;
σ[3] = σ₃;



S₁ = (1 / 2) * σ₁;
S₂ = (1 / 2) * σ₂;
S₃ = (1 / 2) * σ₃;
S = cat(S₁, S₂, S₃, dims=3);



