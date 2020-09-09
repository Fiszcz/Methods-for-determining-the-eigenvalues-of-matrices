include("../utils/verifyEigenvalues.jl")
include("./qrMethod.jl")
using LinearAlgebra
using BenchmarkTools
import Random

# Verification of computeQR
B = [1.0 9.0 3.0; 2.0 4.0 2.0; 22.0 2.0 5.0]

(Q, R) = householderQR(B, 3)

println(Q*R) # should be the same as entry matrix

Random.seed!(100)
A = Array(Symmetric(rand(5, 5) * 100))

(eigenvalues, eigenvectors) = qrMethod(A)

display(eigvals(A))

display(eigenvectors)

println(eigenvalues)

display(@benchmark qrMethod(A))
