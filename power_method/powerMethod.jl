using LinearAlgebra

function powerMethod(matrix::Array{Float64,2})
    eigenvector = ones(size(matrix, 2))

    difference::Float64 = 0.001
    lastEigenValue::Float64 = 0.0
    newEigenValue::Float64 = 0.0

    while true
        eigenvector = normalize(matrix * eigenvector)

        transposed = eigenvector'
        numeratorValue = @fastmath (transposed * matrix * eigenvector)
        denominatorValue = @fastmath (transposed * eigenvector)

        newEigenValue = @fastmath numeratorValue / denominatorValue

        if (abs(newEigenValue - lastEigenValue) <= difference)
            break
        end

        lastEigenValue = newEigenValue
    end

    return (newEigenValue, eigenvector)
end

function powerMethodWithDeflation(matrix::Array{Float64, 2})
    matrixDimension = size(matrix, 2)

    eigenvalues = Vector{Float64}(undef, matrixDimension)
    eigenvectors = Array{Float64, 2}(undef, matrixDimension, matrixDimension)
    for i in eachindex(eigenvalues)
        (newEigenvalue, newEigenvector) = powerMethod(matrix)
        @inbounds eigenvalues[i] = newEigenvalue
        @inbounds eigenvectors[:, i] = newEigenvector
        if i != matrixDimension
            @fastmath matrix -= newEigenvalue*newEigenvector*newEigenvector'
        end
    end
    return (eigenvalues, eigenvectors)
end
