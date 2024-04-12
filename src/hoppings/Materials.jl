⊗(A, B) = kron(A, B)
×(u, v) = cross(u, v)

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


mutable struct Hopping
    a::Int # orbital/site index 1
    b::Int # orbital/site index 2 with PBC
    ia # index vector of site A
    ib # index vector of site B without PBC
    ra # location of atom a
    rb # location of atom b
    r # radius from a to b
    t  # hopping parameter affiliated with c†₂c₁ in spin basis. (i.e., t*I(2) or t*σ₊ might make sense)
    edge::Bool # does this hop off the edge of the superlattice?
    N # vector describing the [n₁;n₂;n₃]⋅[a₁;a₂;a₃] superlattice unit cell of site ib
    desc::String
end

mutable struct Electrode
    xrange::Vector{Int}
    yrange::Vector{Int}
    zrange::Vector{Int}
    n::Int
    connectfrom::String # direction connecting from (i.e, "-x","+x",etc)
    type::String # material that the electrode is made of
    A::Function # magnitude of the exchange field
end

function metalHopping(p::Dict, NNs::Vector{Hopping}, ia::Vector{Int})
    iorb = ia[5]
    t = (3 * p["t"] - 1 * eV) * (I(2))
    pushHopping!(NNs, t, ia, ia, p)
    for ax = 1:3
        for dir = [-1, 1]
            # for weyl term in hamiltonian
            di = zeros(5)
            di[ax] = dir
            ib = Int.(ia + di)
            #Ra = xyztor(p,ia); Rb = xyztor(p,ib); 
            #δ = Rb - Ra
            # implement H = +vf*𝐩⋅𝛔 = -vf𝑖ħ ∇ᵣ⋅σ on finite grid
            t = -1 * p["t"] * I(2)
            pushHopping!(NNs, t, ia, ib, p)
            #t = (p.ϵ₁ + 2*p.t)*(I(2))
            #pushHopping!(NNs, t, ia, ia, p)
        end
    end
end

function insHopping(p::Dict, NNs::Vector{Hopping}, ia::Vector{Int})
    iorb = ia[5]
    t = (p["ε₁"] + 3 * p["t"]) * (I(2))
    pushHopping!(NNs, t, ia, ia, p)
    for ax = 1:3
        for dir = [-1, 1]
            # for weyl term in hamiltonian
            di = zeros(5)
            di[ax] = dir
            ib = Int.(ia + di)
            t = -1 / 2 * p["t"] * I(2)
            pushHopping!(NNs, t, ia, ib, p)
            #t = (p.ϵ₁ + 2*p.t)*(I(2))
            #pushHopping!(NNs, t, ia, ia, p)
        end
    end
end

function nextsite(iorb::Int)
    return (1 - 2 * iorb)
end


function pushHopping!(NNs::Vector, t, ia::Vector{Int}, ib::Vector{Int}, p)
    a = xyztoi(p, ia)
    b = xyztoi(p, ib)
    ra = xyztor(p, ia)
    rb = xyztor(p, ib)
    r = rb - ra
    # for hopping term
    NN = deepcopy(Hopping(a, b, ia, ib, ra, rb, r, t, false, [0; 0; 0], ""))
    push!(NNs, NN)
end

function weyl3Hopping(p::Dict, NNs::Vector{Hopping}, ia::Vector{Int})
    iorb = ia[5]
    ib = ia
    # nearest neighbors
    for ax = 1:3
        for dir = [-1, 1]
            # for weyl term in hamiltonian
            di = zeros(5)
            di[ax] = dir
            #di[5] = nextsite(iorb); 
            ib = Int.(ia + di)
            t = zeros(ComplexF64, 2, 2)
            if (ax == 3)
                t = p["t"] * σ[3]
            elseif (ax == 1)
                t == p["t"] * (-σ[3] .+ dir * 2 * im * σ[1])
            elseif (ax == 2)
                t == p["t"] * (-σ[3] .+ dir * 2 * im * σ[2])
            else
                println("Something broken! No 4th axis")
            end
            pushHopping!(NNs, t, ia, ib, p)
        end
    end
    # next nearest neighbors
    for dx = [-1, 1]
        for dy = [-1, 1]
            di = zeros(5)
            di[1] = dx
            di[2] = dy
            ib = Int.(ia + di)
            t = 1.5 * im * p["t"] * (-dx * σ[1] .+ -dy * σ[2]) # implements the next-nearest neigbor term
            pushHopping!(NNs, t, ia, ib, p)
        end
    end
    # next-to-next nearest neighbors
    for ax = 1:2
        for dir = [-1, 1]
            di = zeros(5)
            di[ax] = 2 * dir
            ib = Int.(ia + di)
            t = dir * p["t"] * 0.5 * im * σ[ax]
            pushHopping!(NNs, t, ia, ib, p)
        end
    end
end


function weyl2Hopping(p::Dict, NNs::Vector{Hopping}, ia::Vector{Int})
    iorb = ia[5]
    ib = ia
    # nearest neighbors
    for ax = 1:3
        for dir = [-1, 1]
            # for weyl term in hamiltonian
            di = zeros(5)
            di[ax] = dir
            #di[5] = nextsite(iorb); 
            ib = Int.(ia + di)
            t = zeros(ComplexF64, 2, 2)
            if (ax == 1)
                t = p["t"] * (σ[1] .- σ[3])
            elseif (ax == 2)
                t = -p["t"] * (σ[1] .+ σ[3])
            elseif (ax == 3)
                t = p["t"] * σ[3]
            else
                println("Something broken! No 4th axis")
            end
            #t = (-im/2)*dir*p.t*σ[ax]
            pushHopping!(NNs, t, ia, ib, p)
        end
    end
    for dx = [-1, 1]
        for dy = [-1, 1]
            di = zeros(5)
            di[1] = dx
            di[2] = dy
            ib = Int.(ia + di)
            t = -dx * dy * 0.5 * p["t"] * σ[2] # implements the next-nearest neigbor term
            pushHopping!(NNs, t, ia, ib, p)
        end
    end
end

function chern2DHopping(p::Dict, NNs::Vector{Hopping}, ia::Vector{Int})
    iorb = ia[5]
    ib = ia
    #ib[5] += nextsite(iorb)
    #t = 3*p.t*(I(2))
    #pushHopping!(NNs, t, ia, ib , p)
    # for H₂ = τ₃⊗σ₀
    t = nextsite(iorb) * (2 * p["m2"] + p["γ"]) * (I(2))
    pushHopping!(NNs, t, ia, ia, p)
    for ax = 1:3
        for dir = [-1, 1]
            # for weyl term in hamiltonian
            di = zeros(5)
            di[ax] = dir
            # for Hweyl = τ₁⊗k⋅σ
            di[5] = nextsite(iorb)
            ib = Int.(ia + di)
            t = (im / 2) * dir * p["t"] * σ[ax]
            #Ra = xyztor(p,ia); Rb = xyztor(p,ib); 
            #δ = Rb - Ra
            # implement H = +vf*𝐩⋅𝛔 = -vf𝑖ħ ∇ᵣ⋅σ on finite grid
            #t = nextsite(iorb)*(-im/2)*dir*p.t*σ[ax]
            pushHopping!(NNs, t, ia, ib, p)
            # for normal hopping term in hamiltonian
            ib[5] = iorb

            t = -(1 / 2) * nextsite(iorb) * p["m2"] * (I(2))
            pushHopping!(NNs, t, ia, ib, p)
        end
    end
end


function weylHopping(p::Dict, NNs::Vector{Hopping}, ia::Vector{Int})
    iorb = ia[5]
    ib = ia
    #ib[5] += nextsite(iorb)
    #t = 3*p.t*(I(2))
    #pushHopping!(NNs, t, ia, ib , p)
    # for H₂ = τ₃⊗σ₀
    t = 3 * nextsite(iorb) * p["t"] * (I(2))
    pushHopping!(NNs, t, ia, ia, p)
    for ax = 1:3
        for dir = [-1, 1]
            # for weyl term in hamiltonian
            di = zeros(5)
            di[ax] = dir
            # for Hweyl = τ₁⊗k⋅σ
            di[5] = nextsite(iorb)
            ib = Int.(ia + di)
            t = (-im / 2) * dir * p["t"] * σ[ax]
            #Ra = xyztor(p,ia); Rb = xyztor(p,ib); 
            #δ = Rb - Ra
            # implement H = +vf*𝐩⋅𝛔 = -vf𝑖ħ ∇ᵣ⋅σ on finite grid
            #t = nextsite(iorb)*(-im/2)*dir*p.t*σ[ax]
            pushHopping!(NNs, t, ia, ib, p)
            # for normal hopping term in hamiltonian
            ib[5] = iorb

            t = -(1 / 2) * nextsite(iorb) * p["t"] * (I(2))
            pushHopping!(NNs, t, ia, ib, p)
        end
    end
end


hoppingDict = Dict{String,Function}(
    "mtjweyl" => weylHopping,
    "2Dchern" => chern2DHopping,
    "weyl" => weylHopping,
    "weyl2" => weyl2Hopping,
    "weyl3" => weyl3Hopping,
    "wins" => insHopping,
    "ins" => insHopping,
    "insulator" => insHopping,
    "metal" => metalHopping
)


#end