# include( ("../src/UsefulFunctions.jl" ))
using SparseArrays
using SuiteSparse
using LinearAlgebra
using Plots
exmat(n) = sprand(n,n,3/n)+5I(n);
#f(n) = recInv(sprand(n,n,3/n)+5I)
Φ = 0.2001
# Φ = 0.000

# Create the sparse matrix
H(n) = sparse(Tridiagonal(exp(im*Φ)*ones(4*n-1), ComplexF64.(3:1:(4*(n)+2)), exp(-im*Φ)*ones(4*n-1)))
# H(n) = sparse(Tridiagonal(exp(im*Φ)*ones(4*n-1),2*ComplexF64.(ones(4*n)),exp(-im*Φ)*ones(4*n-1)))
# H(n) = sparse(Tridiagonal(exp(im*Φ)*ones(4*n-1),2*ComplexF64.(ones(4*n)),exp(-im*Φ)*ones(4*n-1)))

f(n) = PartialRecInvGr(H(n), 4, "transport");
g(n) = RowsGrInv(H(n), 4);



function b(row::Int,B::Int=4)
    return ((row-1)*B+1):(row*B)
end

function buildM!(Mrows::Vector{Int},Mcols::Vector{Int},elems::Vector{ComplexF64}, rows::UnitRange{Int},cols::UnitRange{Int},submat)
    offsetrow = rows[1]; offsetcol = cols[1]
    for row in rows
        for col in cols
            push!(Mrows,row); push!(Mcols,col)
            push!(elems,submat[row-offsetrow+1,col-offsetcol+1])
        end
    end
end    

function labeldisplay(submat, label::String)
    println("$label:")
    display(submat)
end

function exterior(gʳp::Matrix, V::SparseMatrixCSC, Gʳprev::Matrix)
    return transpose(gʳp*V*transpose(Gʳprev))
end
   


function PartialRecInvGr1(M::SparseMatrixCSC, B::Int=1, offdiag = true)

    n = size(M)[1];
    nBs = Int(n/B) # number of diagonal blocks
    LCblocks = Matrix[] # go down the diagonal from left = backward generator but testing shows opposite
    RCblocks = Matrix[] # go down the diagonal from right = forward generator

    
    # final inverted matrix
    Minv = spzeros(ComplexF64,n,n)
    rows = Int[]; cols = Int[]; elems = ComplexF64[];
    gʳL₀ = inv(Array(M[1:B,1:B]))
    # println(gʳL₀)


    # diagonal block hamiltonian
    push!(LCblocks,gʳL₀)
    for ib = 2:nBs
        gʳLprior = LCblocks[ib-1]
        # off-diagonal block coupling hamiltonian
        V = M[b(ib-1,B),b(ib,B)]
        #display(V)
        # on-diagonal block hamiltonian
        Dᵢ = M[b(ib,B),b(ib,B)]
        # effective surface green's function including this point
        gʳLᵢ = inv(Array(Dᵢ - V'*gʳLprior*V))
        # println(gʳLᵢ)

        push!(LCblocks,copy(gʳLᵢ))
    end


    # now we will go back up the diagonal and incorporate coupling from right
    # println("LCblocks:", LCblocks)
    Gʳend = last(LCblocks)
    Gʳᵢ = Gʳend
    buildM!(rows,cols,elems,b(nBs,B),b(nBs,B),Gʳend)
    Gʳplus = copy(Gʳend)

    println(LCblocks)


    println(Gʳend)
    for ib = reverse(1:(nBs-1))
        gʳLᵢ = LCblocks[ib]
        # off-diagonal block coupling hamiltonian
        # top off diagonal
        V = M[b(ib,B),b(ib+1,B)]

        # calculating diagonals
        Gʳᵢ = gʳLᵢ*(I(B) + V*Gʳplus*V'*gʳLᵢ)
        println(Gʳᵢ)
        buildM!(rows,cols,elems,b(ib,B),b(ib,B),Array(Gʳᵢ))

        #Gʳplus = last diagonal block 
        #gʳLᵢ = forward/backward generator sequence
        Gʳplus = copy(Gʳᵢ)

    end
    return sparse(rows,cols,elems)

end

function PartialRecInvGr(M::SparseMatrixCSC, B::Int=1, offdiag = true)
    n = size(M)[1];
    nBs = Int(n/B) # number of diagonal blocks
    LCblocks = Matrix[] # go down the diagonal from left = backward generator but testing shows opposite
    RCblocks = Matrix[] # go down the diagonal from right = forward generator

    
    # final inverted matrix
    Minv = spzeros(ComplexF64,n,n)
    rows = Int[]; cols = Int[]; elems = ComplexF64[];
    gʳL₀ = inv(Array(M[1:B,1:B]))
    println(gʳL₀)


    # diagonal block hamiltonian
    push!(LCblocks,gʳL₀)
    for ib = 2:nBs
        gʳLprior = LCblocks[ib-1]
        # off-diagonal block coupling hamiltonian
        V = M[b(ib-1,B),b(ib,B)]
        #display(V)
        # on-diagonal block hamiltonian
        Dᵢ = M[b(ib,B),b(ib,B)]
        # effective surface green's function including this point
        gʳLᵢ = inv(Array(Dᵢ - V'*gʳLprior*V))
        println(gʳLᵢ)

        push!(LCblocks,copy(gʳLᵢ))
    end


    # now we will go back up the diagonal and incorporate coupling from right
    # println("LCblocks:", LCblocks)
    Gʳend = last(LCblocks)
    Gʳᵢ = Gʳend
    buildM!(rows,cols,elems,b(nBs,B),b(nBs,B),Gʳend)
    Gʳplus = copy(Gʳend)



    Gʳbotplus = copy(Gʳend)
    for ib = reverse(1:(nBs-1))
        gʳLᵢ = LCblocks[ib]
        # off-diagonal block coupling hamiltonian
        V = M[b(ib,B),b(ib+1,B)]
        # see Klimeck's pwpt
        # now do the block-diagonal Gʳ
        Gʳᵢ = gʳLᵢ*(I(B) + V*Gʳplus*V'*gʳLᵢ)
        #Atruth(ib,ib,Gʳᵢ,"diagonal subblock")
        buildM!(rows,cols,elems,b(ib,B),b(ib,B),Array(Gʳᵢ));
        Gʳᵢbot = exterior(gʳLᵢ,V,Gʳbotplus)
        #truth(nBs,ib,Gʳᵢbot,"offdiagonal block")
        #Gʳᵢbot = transpose(-gʳLᵢ*conj.(V)*transpose(Gʳbotplus))
        buildM!(rows,cols,elems,b(nBs,B),b(ib,B),Gʳᵢbot)
        Gʳplus = copy(Gʳᵢ)
        Gʳbotplus = copy(Gʳᵢbot)
    end
    # generate the right-connected gʳs
    gʳR₀ = inv(Array(M[b(nBs,B),b(nBs,B)]))
    push!(RCblocks,copy(gʳR₀))
    gʳRprior = gʳR₀
    for ib = reverse(1:(nBs-1))
        # off-diagonal block coupling hamiltonian
        V = M[b(ib,B),b(ib+1,B)]
        #display(V)+
        # on-diagonal block hamiltonian
        Dᵢ = M[b(ib,B),b(ib,B)]
        # effective surface green's function including this point
        gʳRᵢ = inv(Array(Dᵢ - V*gʳRprior*V'))
        push!(RCblocks,copy(gʳRᵢ))
        gʳRprior = gʳRᵢ
    end

    # and now do the top row
    reverse!(RCblocks)

    println("LCblocks:", LCblocks)
    println("RCblocks:", RCblocks)

    Gʳtopmin = Gʳᵢ
    println("Gʳᵢ", Gʳᵢ)
    for ib = 2:nBs
        gʳRᵢ = RCblocks[ib]
        V = M[b(ib,B),b(ib-1,B)]
        Gʳᵢtop = exterior(gʳRᵢ,V,Gʳtopmin)
        #Gʳᵢtop = -transpose(gʳRᵢ*conj.(V)*transpose(Gʳtopmin))
        buildM!(rows,cols,elems,b(1,B),b(ib,B),Gʳᵢtop)
        Gʳtopmin = copy(Gʳᵢtop)
    end
    
    
    return sparse(rows,cols,elems)
end


dGʳ(n) = PartialRecInvGr(H(n), 1, "DOS")
fullGʳ(n) = inv(Array(H(n)))

# println(Array(H(2)))

aa = fullGʳ(2)
# bb = Matrix(dGʳ(2))
# println(size(aa))

# for i in 1:8
#     print(aa[i, i])
#     print("    ")
#     print(bb[i, i])
#     println()
# end
