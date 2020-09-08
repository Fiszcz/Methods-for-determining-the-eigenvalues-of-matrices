using LinearAlgebra

function powerMethod(matrix::Array{Float64, 2})
    eigenVector = ones(size(matrix, 2))

    difference::Float64 = 0.001
    lastEigenValue::Float64 = 0
    newEigenValue::Float64 = 0

    while true
        lastEigenValue = newEigenValue

        eigenVector = normalize(matrix * eigenVector)

        transposed = eigenVector'
        numeratorValue = (transposed * matrix * eigenVector)
        denominatorValue = (transposed * eigenVector)

        newEigenValue = numeratorValue / denominatorValue

        if (abs(newEigenValue - lastEigenValue) <= difference)
            break
        end
    end

    return (newEigenValue, eigenVector)
end

function powerMethodWithDeflation(matrix::Array{Float64, 2})
    matrixDimension = size(matrix, 2)
    eigenvalues = Vector{Float64}()
    eigenvectors = []
    for i = 1:matrixDimension
        (newEigenvalue, newEigenvector) = powerMethod(matrix)
        append!(eigenvalues, newEigenvalue)
        append!(eigenvectors, [newEigenvector])
        if i != matrixDimension
            matrix = matrix - newEigenvalue*newEigenvector*newEigenvector'
        end
    end
    return (eigenvalues, eigenvectors)
end

A = [2.0 1.0; 1.0 2.0]
B = [5.0 3.0; 3.0 5.0]
C = [1.0 2.0 3.0; 2.0 4.0 2.0; 3.0 2.0 5.0]

powerMethod(C)

display(eigen(C))
display(powerMethodWithDeflation(C))
