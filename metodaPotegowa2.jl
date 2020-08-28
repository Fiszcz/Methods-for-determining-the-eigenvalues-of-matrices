using LinearAlgebra

function powerIteration(matrix::Array{Float64, 2}, size::Int64)
    eigenVector = ones(size)

    difference::Float64 = 0.001
    lastEigenValue::Float64 = 0
    newEigenValue::Float64 = 0

    while true
        lastEigenValue = newEigenValue

        eigenVector = normalize(matrix * eigenVector)

        transposed = eigenVector'
        avalorNumerador = (transposed * matrix * eigenVector)
        avalorDenominador = (transposed * eigenVector)

        newEigenValue = avalorNumerador / avalorDenominador

        if (abs(newEigenValue - lastEigenValue) <= difference)
            break
        end
    end

    println(newEigenValue)
    println(eigenVector)

end

A = [2.0 1.0; 1.0 2.0]
B = [5.0 3.0; 3.0 5.0]
C = [1.0 2.0 3.0; 2.0 4.0 2.0; 3.0 2.0 5.0]

powerIteration(C, 3)
