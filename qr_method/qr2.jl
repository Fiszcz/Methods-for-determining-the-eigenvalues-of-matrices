include("./verifyEigenvalues.jl")

using LinearAlgebra

function computeQR(matrix::Matrix{Float64}, matrixDimension::Int64)
    R = zeros(Float64, matrixDimension, matrixDimension)
    Q = copy(R)

    for j = 1:matrixDimension
        v = matrix[:,j]
        for i = 1:j-1
            R[i, j] = Q[:,i]'*v
            v = v - R[i,j]*Q[:,i]
        end
        R[j,j] = norm(v)
        Q[:,j] = v / R[j,j]
    end

    return (Q, R)
end
# Verification of computeQR
# A = [1.0 9.0 3.0; 2.0 4.0 2.0; 22.0 2.0 5.0]
# (Q, R) = computeQR(A, 3)
# println(Q*R) # should be the same as entry matrix

function getEigenvalues(matrix::Matrix{Float64})
    matrixDimension = size(matrix, 2)

    An = matrix
    α = 0

    identityMatrix = Matrix{Float64}(I, matrixDimension, matrixDimension)

    P = 1

    for i = 1:100
        (Q, R) = computeQR(An - identityMatrix * α, matrixDimension)
        P = P * Q
        An = R * Q + identityMatrix * α
    end

    return (diag(An), P)
end

A = [2.0 3.0 5.0; 3.0 8.0 1.0; 5.0 1.0 3.0]

(eigenvalues, eigenvectors) = getEigenvalues(A)

println(eigenvalues)

verifyEigenvalues(A, eigenvalues)

println([eigenvectors[:, i] for i in 1:size(A, 2)])

verifyEigenvectors(A, [eigenvectors[:, i] for i in 1:size(A, 2)])
