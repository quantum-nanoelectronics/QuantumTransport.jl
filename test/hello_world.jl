# this file includes a test for a sample Hello World test

module TestHelloWorld
using QuantumTransport
using Test

println("\033[1mRunning Hello World Test\033[0m")

@test sampleTest() == "Hello World!"

end