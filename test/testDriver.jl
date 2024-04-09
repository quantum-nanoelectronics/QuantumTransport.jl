module testDriver
using QuantumTransport
include("InBi.jl")
include("../src/driver/Driver.jl")


function DriverTest(pNNs::NamedTuple, A::Function)
    dict = Dict()
    
    for (key,value) in zip(keys(pNNs), pNNs)
        dict[string(key)] = value
    end
    #print(dict)

    main(dict, A)
end

p, p1, p2, p3, A = InBi.generateParams()
DriverTest(p, A)

end