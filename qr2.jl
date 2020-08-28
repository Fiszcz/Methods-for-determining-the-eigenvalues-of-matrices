module QR_Computation

using LinearAlgebra

include("./verifyEigenvalues.jl")

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

function getEigenvalues(matrix::Matrix{Float64})
    matrixDimension = size(matrix, 2)

    An = matrix
    α = 0

    identityMatrix = Matrix{Float64}(I, matrixDimension, matrixDimension)

    for i = 1:100
        (Q, R) = computeQR(An - identityMatrix, matrixDimension)
        An = R * Q + identityMatrix * α
    end

    return diag(An)
end

A = [1.0 9.0 3.0; 2.0 4.0 2.0; 22.0 2.0 5.0]

result = getEigenvalues(A)

println(result)

Main.verifying.verifyEigenvalues(A, result)

end
