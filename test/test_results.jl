include("../power_method/powerMethod.jl")
include("../qr_method/qrMethod.jl")
using LinearAlgebra

function eigenvalues_from_LinearAlgebra(entryMatrix)
    printstyled("Computed values by LinearAlgebra package:\n", color = :(green))
    (values, vectors) = eigen(entryMatrix)
    println("Eigenvalues: ", values)
    println("Eigenvectors: ")
    display(vectors)
end

function test_results_of_power_method(entryMatrix)
    (eigenValues, eigenVectors) = powerMethodWithDeflation(entryMatrix)

    printstyled("POWER METHOD WITH DEFLATION:\n", color = :(green))

    println("Computed values:")
    println("Eigenvalue: ", eigenValues)
    println("Eigenvector: ")
    display(eigenVectors)
end

function test_results_of_qr_method(entryMatrix)
    (eigenValues, eigenVectors) = qrMethod(entryMatrix)

    printstyled("QR METHOD\n", color = :(green))

    println("Computed values:")
    println("Eigenvalues: ", eigenValues)
    println("Eigenvectors: ")
    display(eigenVectors)
end

entryMatrix = Array(Symmetric(rand(3, 3)))

eigenvalues_from_LinearAlgebra(entryMatrix)
test_results_of_power_method(entryMatrix)
test_results_of_qr_method(entryMatrix)
