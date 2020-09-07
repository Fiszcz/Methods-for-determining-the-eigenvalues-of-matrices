include("../utils/verifyEigenvalues.jl")
using LinearAlgebra

# Verification of computeQR
A = [1.0 9.0 3.0; 2.0 4.0 2.0; 22.0 2.0 5.0]

(Q, R) = computeQR(A, 3)

println(Q*R) # should be the same as entry matrix

B = [2.0 3.0 5.0; 3.0 8.0 1.0; 5.0 1.0 3.0]

(eigenvalues, eigenvectors) = getEigenvalues(B)

println(eigenvalues)

verifyEigenvalues(B, eigenvalues)

println([eigenvectors[:, i] for i in 1:size(B, 2)])

verifyEigenvectors(B, [eigenvectors[:, i] for i in 1:size(B, 2)])
