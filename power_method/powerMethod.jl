using LinearAlgebra

function powerMethod(matrix::Symmetric{Float64,Array{Float64,2}})
    eigenVector = ones(size(matrix, 2))

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

    return (newEigenValue, eigenVector)
end
