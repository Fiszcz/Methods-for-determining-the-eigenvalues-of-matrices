include("../jacobi_method/jacobiMethod.jl")
include("../power_method/powerMethod.jl")
include("../qr_method/qrMethod.jl")
using LinearAlgebra

function eigenvalues_from_LinearAlgebra(entryMatrix)
    printstyled("Computed values by LinearAlgebra package:\n", color = :(green))
    (values, vectors) = eigen(entryMatrix)
    println("Eigenvalues: ", values)
    println("Eigenvectors: ")
    display(vectors)
    println()
    println()
end

function test_results_of_jacobi_method(entryMatrix)
    (eigenValue, eigenVector) = jacobiMethod(entryMatrix)

    printstyled("JACOBI METHOD:\n", color = :(green))

    println("Computed values:")
    println("Eigenvalue: ", eigenValue)
    println("Eigenvector: ")
    display(eigenVector)
    println()
    println()
end

function test_results_of_power_method(entryMatrix)
    (eigenValue, eigenVector) = powerMethod(entryMatrix)

    printstyled("POWER METHOD:\n", color = :(green))

    println("Computed values:")
    println("Eigenvalue: ", eigenValue)
    println("Eigenvector: ")
    display(eigenVector)
    println()
end

function test_results_of_qr_method(entryMatrix)
    (eigenValues, eigenVectors) = qrMethod(entryMatrix)

    printstyled("QR METHOD\n", color = :(green))

    println("Computed values:")
    println("Eigenvalues: ", eigenValues)
    println("Eigenvectors: ")
    display(eigenVectors)
    println()
    println()
end

entryMatrix = Symmetric(rand(5, 5))

eigenvalues_from_LinearAlgebra(entryMatrix)
test_results_of_jacobi_method(entryMatrix)
test_results_of_power_method(entryMatrix)
test_results_of_qr_method(entryMatrix)
