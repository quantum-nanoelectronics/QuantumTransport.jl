"""
This is the main test file that is called first when running any tests.
Uses modules in the test files to avoid polluting the global namespace, with QuantumTransport, Test, and other packages used. 
"""

include("runtests_unit.jl")
include("run_qt.jl")
