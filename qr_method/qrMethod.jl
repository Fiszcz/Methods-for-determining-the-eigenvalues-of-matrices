using LinearAlgebra

function compute_Q(W::Array{Float64,2}, matrixDimension::Int64)
    Z = Array{Float64,2}(I, matrixDimension, matrixDimension)
    @views for column = matrixDimension:-1:1
        Z -= W[:,column]*(W[:,column]'*Z)
    end
    return Z
end

const sqrt_2 = sqrt(2)
function householder_reflection(x)
    nw = norm(x)
    if nw != 0
        w = x/nw
        w[1] += (w[1] == 0.0 ? 1.0 : sign(w[1]))
        w /= sqrt(abs(w[1]))
    else
        w = x
        w[1] = sqrt_2
    end

    return w
end

function householderQR(matrix::Array{Float64,2}, matrixDimension::Int64)
    W = zeros(Float64, matrixDimension, matrixDimension)
    R::Array{Float64,2} = copy(matrix)

    for j = 1:matrixDimension
        w = householder_reflection(@view(R[j:end, j]))
        W[j:end, j] = w
        R[j:end, j:end] -= w*(w'*@view(R[j:end, j:end]))
        R[j+1:end, j] .= 0.0
    end
    Q::Array{Float64, 2} = compute_Q(W, matrixDimension)

    return (Q, R)
end

function qrMethod(matrix::Array{Float64,2})
    matrixDimension = size(matrix, 2)

    α = matrix[end, end]

    identityMatrix = Array{Float64,2}(I, matrixDimension, matrixDimension)

    P = Array{Float64,2}(I, matrixDimension, matrixDimension)

    previous_diagonal_length::Float64 = 10.0
    current_diagonal_length_squared::Float64 = eps(Float64)
    while abs((previous_diagonal_length - current_diagonal_length_squared)/previous_diagonal_length) > 0.001
        (Q, R) = householderQR(matrix - identityMatrix * α, matrixDimension)
        P *= Q
        matrix = R * Q + identityMatrix * α

        previous_diagonal_length = current_diagonal_length_squared
        current_diagonal_length_squared = sum(diag(matrix).^2)
    end

    return (diag(matrix), P)
end
