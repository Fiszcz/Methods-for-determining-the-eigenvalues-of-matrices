using LinearAlgebra

function computeQR(matrix::Array{Float64,2}, matrixDimension::Int64)
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

function qrMethod(matrix::Symmetric{Float64,Array{Float64,2}})
    matrixDimension = size(matrix, 2)

    An = matrix
    α = 0

    identityMatrix = Array{Float64,2}(I, matrixDimension, matrixDimension)

    P = 1

    for i = 1:100
        (Q, R) = computeQR(An - identityMatrix * α, matrixDimension)
        P = P * Q
        An = R * Q + identityMatrix * α
    end

    return (diag(An), P)
end
