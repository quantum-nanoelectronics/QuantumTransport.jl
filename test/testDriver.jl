module testDriver
using QuantumTransport
using Test

include("InBi.jl")

# TODO Remove include statements from /src folder in test files like this one (except do it for all test files)
include("../src/driver/Driver.jl")


function DriverTest(pNNs::NamedTuple, A::Function)
    dict = Dict()
    
    for (key,value) in zip(keys(pNNs), pNNs)
        dict[string(key)] = value
    end
    # println(dict)

    main(dict, A)

    return true
end

# ϵ₁
# ε₁


p, p1, p2, p3, A = InBi.generateParams()

# added ε₁ = 2.0,
# anotherP = (t = 1.0, t₁ = 0.1, t₂ = 0.4, t₃ = 0.2, t₄ = 0, t₅ = 0.0, t₆ = 0.0, t₇ = 0.0, t₈ = 0.1, t₉ = 0.0, vf = 1000000, η = 0.0001, ε = 0.0, ε₁ = 2.0, ϵ₁ = 2.0, a₁ = [1.0e-9, 0.0, 0.0], a₂ = [0.0, 1.0e-9, 0.0], a₃ = [0.0, 0.0, 1.0e-9], A = [1.0000000000000001e-7 0.0 0.0; 0.0 1.0e-9 0.0; 0.0 0.0 1.0e-9], a = 1.0e-9, b = 1.0e-9, c = 1.0e-9, SLa₁ = [1.0000000000000001e-7, 0.0, 0.0], SLa₂ = [0.0, 1.0e-9, 0.0], SLa₃ = [0.0, 0.0, 1.0e-9], nx = 100, ny = 1, nz = 1, n = 100, norb = 2, nsite = 1, kdict = Dict{String, Vector}("Z" => [0.0, 0.0, 3.1415926535897927e9], "-X₂" => [0.0, -3.1415926535897927e9, 0.0], "-X₃" => [0.0, 0.0, -3.1415926535897927e9], "Γ + 10*iX₁" => ComplexF64[0.0 + 3.1415926535897934e8im, 0.0 + 0.0im, 0.0 + 0.0im], "M" => [3.1415926535897933e7, 3.1415926535897927e9, 0.0], "-Z" => [0.0, 0.0, -3.1415926535897927e9], "X₂" => [0.0, 3.1415926535897927e9, 0.0], "A" => [3.1415926535897933e7, 3.1415926535897927e9, 3.1415926535897927e9], "X₁" => [3.1415926535897933e7, 0.0, 0.0], "Γ + iX₁" => ComplexF64[0.0 + 3.1415926535897933e7im, 0.0 + 0.0im, 0.0 + 0.0im], "-X₁" => [-3.1415926535897933e7, 0.0, 0.0], "X₃" => [0.0, 0.0, 3.1415926535897927e9], "Γ" => [0.0, 0.0, 0.0]), μ = 0.0, μ_disorder = 0.0, E_samples = [0.1], nk = 0, δV = 0.01, T = 300.0, path = "../outputs/testrunstacc/Wed-10-Apr-2024--07.45.20/bloch/", savedata = true, save = true, β = 0.25, runtype = "multiblochdws", fieldtype = "β", ηD = 0.0001, l_scattering = 0.0, parallel = "k", n_BLAS = 8, transport = true, verbose = false, plotfield = true, bands = false, mixedDOS = false, θ = 360.0, sweep = "plot", returnvals = ["transmission"], electrodeMagnetization = true, electrodeMaterial = "weyl", deviceMagnetization = true, deviceMaterial = "weyl", startDWs = 4.0e-8, DWwidth = 9.000000000000001e-9, DWspacing = 1.5000000000000002e-8, λ = 1.0e7, prune = Any[], B = [6.2831853071795866e7 0.0 0.0; 0.0 6.283185307179585e9 0.0; 0.0 0.0 6.283185307179585e9], arpack = true, klist = ["M", "Γ", "X₁", "M", "X₂", "Γ", "X₃"])
# anotherP1 = (t = 1.0, t₁ = 0.1, t₂ = 0.4, t₃ = 0.2, t₄ = 0, t₅ = 0.0, t₆ = 0.0, t₇ = 0.0, t₈ = 0.1, t₉ = 0.0, vf = 1000000, η = 0.0001, ε = 0.0, ε₁ = 2.0, ϵ₁ = 2.0, a₁ = [1.0e-9, 0.0, 0.0], a₂ = [0.0, 1.0e-9, 0.0], a₃ = [0.0, 0.0, 1.0e-9], A = [1.0000000000000001e-7 0.0 0.0; 0.0 1.0e-9 0.0; 0.0 0.0 1.0e-9], a = 1.0e-9, b = 1.0e-9, c = 1.0e-9, SLa₁ = [1.0000000000000001e-7, 0.0, 0.0], SLa₂ = [0.0, 1.0e-9, 0.0], SLa₃ = [0.0, 0.0, 1.0e-9], nx = 100, ny = 1, nz = 1, n = 100, norb = 2, nsite = 1, kdict = Dict{String, Vector}("Z" => [0.0, 0.0, 3.1415926535897927e9], "-X₂" => [0.0, -3.1415926535897927e9, 0.0], "-X₃" => [0.0, 0.0, -3.1415926535897927e9], "Γ + 10*iX₁" => ComplexF64[0.0 + 3.1415926535897934e8im, 0.0 + 0.0im, 0.0 + 0.0im], "M" => [3.1415926535897933e7, 3.1415926535897927e9, 0.0], "-Z" => [0.0, 0.0, -3.1415926535897927e9], "X₂" => [0.0, 3.1415926535897927e9, 0.0], "A" => [3.1415926535897933e7, 3.1415926535897927e9, 3.1415926535897927e9], "X₁" => [3.1415926535897933e7, 0.0, 0.0], "Γ + iX₁" => ComplexF64[0.0 + 3.1415926535897933e7im, 0.0 + 0.0im, 0.0 + 0.0im], "-X₁" => [-3.1415926535897933e7, 0.0, 0.0], "X₃" => [0.0, 0.0, 3.1415926535897927e9], "Γ" => [0.0, 0.0, 0.0]), μ = 0.0, μ_disorder = 0.0, E_samples = [0.1], nk = 0, δV = 0.01, T = 300.0, path = "../outputs/testrunstacc/Wed-10-Apr-2024--09.41.13/bloch/", savedata = true, save = true, β = 0.25, runtype = "multiblochdws", fieldtype = "β", ηD = 0.0001, l_scattering = 0.0, parallel = "k", n_BLAS = 8, transport = true, verbose = false, plotfield = true, bands = false, mixedDOS = false, θ = 360.0, sweep = "none", returnvals = ["transmission"], electrodeMagnetization = true, electrodeMaterial = "weyl", deviceMagnetization = true, deviceMaterial = "weyl", startDWs = 4.0e-8, DWwidth = 9.000000000000001e-9, DWspacing = 1.5000000000000002e-8, λ = 1.0e7, prune = Any["x"], B = [6.2831853071795866e7 0.0 0.0; 0.0 6.283185307179585e9 0.0; 0.0 0.0 6.283185307179585e9], arpack = true, klist = ["M", "Γ", "X₁", "M", "X₂", "Γ", "X₃"], nelectrodes = 2)
# anotherP3 = (t = 1.0, t₁ = 0.1, t₂ = 0.4, t₃ = 0.2, t₄ = 0, t₅ = 0.0, t₆ = 0.0, t₇ = 0.0, t₈ = 0.1, t₉ = 0.0, vf = 1000000, η = 0.0001, ε = 0.0, ε₁ = 2.0, ϵ₁ = 2.0, a₁ = [1.0e-9, 0.0, 0.0], a₂ = [0.0, 1.0e-9, 0.0], a₃ = [0.0, 0.0, 1.0e-9], A = [1.0000000000000001e-7 0.0 0.0; 0.0 1.0e-9 0.0; 0.0 0.0 1.0e-9], a = 1.0e-9, b = 1.0e-9, c = 1.0e-9, SLa₁ = [1.0000000000000001e-7, 0.0, 0.0], SLa₂ = [0.0, 1.0e-9, 0.0], SLa₃ = [0.0, 0.0, 1.0e-9], nx = 100, ny = 1, nz = 1, n = 100, norb = 2, nsite = 1, kdict = Dict{String, Vector}("Z" => [0.0, 0.0, 3.1415926535897927e9], "-X₂" => [0.0, -3.1415926535897927e9, 0.0], "-X₃" => [0.0, 0.0, -3.1415926535897927e9], "Γ + 10*iX₁" => ComplexF64[0.0 + 3.1415926535897934e8im, 0.0 + 0.0im, 0.0 + 0.0im], "M" => [3.1415926535897933e7, 3.1415926535897927e9, 0.0], "-Z" => [0.0, 0.0, -3.1415926535897927e9], "X₂" => [0.0, 3.1415926535897927e9, 0.0], "A" => [3.1415926535897933e7, 3.1415926535897927e9, 3.1415926535897927e9], "X₁" => [3.1415926535897933e7, 0.0, 0.0], "Γ + iX₁" => ComplexF64[0.0 + 3.1415926535897933e7im, 0.0 + 0.0im, 0.0 + 0.0im], "-X₁" => [-3.1415926535897933e7, 0.0, 0.0], "X₃" => [0.0, 0.0, 3.1415926535897927e9], "Γ" => [0.0, 0.0, 0.0]), μ = 0.0, μ_disorder = 0.0, E_samples = [0.1], nk = 0, δV = 0.01, T = 300.0, path = "../outputs/testrunstacc/Wed-10-Apr-2024--21.05.40/neel/", savedata = true, save = true, β = 0.25, runtype = "multineeldws", fieldtype = "β", ηD = 0.0001, l_scattering = 0.0, parallel = "k", n_BLAS = 8, transport = true, verbose = false, plotfield = true, bands = false, mixedDOS = false, θ = 360.0, sweep = "plot", returnvals = ["transmission"], electrodeMagnetization = true, electrodeMaterial = "weyl", deviceMagnetization = true, deviceMaterial = "weyl", startDWs = 7.200000000000001e-8, DWwidth = 9.000000000000001e-9, DWspacing = 1.5000000000000002e-8, λ = 1.0e7, prune = Any[], B = [6.2831853071795866e7 0.0 0.0; 0.0 6.283185307179585e9 0.0; 0.0 0.0 6.283185307179585e9], arpack = true, klist = ["M", "Γ", "X₁", "M", "X₂", "Γ", "X₃"])
# anotherP1 = merge(p, anotherP)

@test DriverTest(p, A)
# @run DriverTest(anotherP1, A)

# @test DriverTest(p, A)

end