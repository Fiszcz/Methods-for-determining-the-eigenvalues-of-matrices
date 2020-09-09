using LinearAlgebra

function powerMethod(matrix::Array{Float64,2})
    eigenVector = ones(size(matrix, 2))

    difference::Float64 = 0.001
    lastEigenValue::Float64 = 0
    newEigenValue::Float64 = 0

    while true
        eigenVector = normalize(matrix * eigenVector)

        transposed = eigenVector'
        numeratorValue = (transposed * matrix * eigenVector)
        denominatorValue = (transposed * eigenVector)

        newEigenValue = numeratorValue / denominatorValue

        if (abs(newEigenValue - lastEigenValue) <= difference)
            break
        end

        lastEigenValue = newEigenValue
    end

    return (newEigenValue, eigenVector)
end

function powerMethodWithDeflation(matrix::Array{Float64, 2})
    matrixDimension = size(matrix, 2)

    eigenvalues = Vector{Float64}(undef, 5)
    eigenvectors = Array{Float64, 2}(undef, matrixDimension, matrixDimension)
    for i = 1:matrixDimension
        (newEigenvalue, newEigenvector) = powerMethod(matrix)
        eigenvalues[i] = newEigenvalue
        eigenvectors[:, i] = newEigenvector
        if i != matrixDimension
            matrix = matrix - newEigenvalue*newEigenvector*newEigenvector'
        end
    end
    return (eigenvalues, eigenvectors)
end
