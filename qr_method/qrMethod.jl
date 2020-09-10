using LinearAlgebra

@inline function compute_Q(W::Array{Float64,2}, matrixDimension::Int64)
    Z = Array{Float64,2}(I, matrixDimension, matrixDimension)
    @views for column = matrixDimension:-1:1
        Z -= W[:,column]*(W[:,column]'*Z)
    end
    return Z
end

const sqrt_2 = sqrt(2)
@inline function householder_reflection(x)
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

@inline function householderQR(matrix::Array{Float64,2}, matrixDimension::Int64)
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

    inProgress = true
    currently_checked_column = 1
    currently_checked_row = 2
    while inProgress
        (Q, R) = householderQR(matrix - identityMatrix * α, matrixDimension)
        P *= Q
        matrix = R * Q + identityMatrix * α

        inProgress = false
        for column = currently_checked_column:matrixDimension-1
            for row = currently_checked_row:matrixDimension
                if abs(matrix[row, column]) > 0.001
                    inProgress = true
                    currently_checked_column = column
                    currently_checked_row = row
                    break
                end
            end
            if inProgress
                break
            end
            currently_checked_row = column + 2
        end
    end

    return (diag(matrix), P)
end
