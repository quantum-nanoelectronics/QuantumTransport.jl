using QuantumTransport
using Test

# this file includes a test that checks for correctness for the block matrix inversion methods

# Test the woodbury inverse method
@test QuantumTransport.block_inv_main()

# Test the RGF inverse method for correctness only
@test QuantumTransport.rgf_main()