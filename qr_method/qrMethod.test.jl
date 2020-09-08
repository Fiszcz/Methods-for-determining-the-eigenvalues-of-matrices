include("../utils/verifyEigenvalues.jl")
include("./qrMethod0.jl")
using LinearAlgebra

# Verification of computeQR
A = [1.0 9.0 3.0; 2.0 4.0 2.0; 22.0 2.0 5.0]

(Q, R) = householderQR(A, 3)

println(Q*R) # should be the same as entry matrix

B = [2.0 3.0 5.0; 3.0 8.0 1.0; 5.0 1.0 3.0]

display(eigen(B))

(eigenvalues, eigenvectors) = qrMethod(B)

display(eigenvectors)

println(eigenvalues)

verifyEigenvalues(B, eigenvalues)
