using LinearAlgebra

function compute_Q(W, matrixDimension)
    Z = I
    for j = matrixDimension:-1:1
        Z = Z - W[:,j]*(W[:,j]'*Z)
    end
    return Z
end

sqrt_2 = sqrt(2)
function householder_reflection(x)
    nw = norm(x)
    if nw != 0
        w = x/nw
        w[1] = w[1] + (w[1] == 0 ? 1 : sign(w[1]))
        w = w/sqrt(abs(w[1]))
    else
        w = x
        w[1] = sqrt_2
    end

    return w
end

function householderQR(matrix::Array{Float64,2}, matrixDimension::Int64)
    W = zeros(Float64, matrixDimension, matrixDimension)
    R = copy(matrix)

    for j = 1:matrixDimension
        w = householder_reflection(R[j:matrixDimension, j])
        W[j:matrixDimension, j] = w
        R[j:matrixDimension, j:matrixDimension] = R[j:matrixDimension, j:matrixDimension] - w*(w'*R[j:matrixDimension, j:matrixDimension])
        R[j+1:matrixDimension, j] .= 0
    end
    Q = compute_Q(W, matrixDimension)

    return (Q, R)
end

function qrMethod(matrix::Array{Float64,2})
    matrixDimension = size(matrix, 2)

    An = matrix
    α = 0

    identityMatrix = Array{Float64,2}(I, matrixDimension, matrixDimension)

    P = 1

    for i = 1:100
        (Q, R) = householderQR(An - identityMatrix * α, matrixDimension)
        P = P * Q
        An = R * Q + identityMatrix * α
    end

    return (diag(An), P)
end
