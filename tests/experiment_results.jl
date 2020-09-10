include("../jacobi_method/jacobiMethod.jl")
include("../power_method/powerMethod.jl")
include("../qr_method/qrMethod.jl")
include("../utils/verifyEigenvalues.jl")
using LinearAlgebra

function eigenvalues_from_LinearAlgebra(entryMatrix)
    printstyled("Computed values by LinearAlgebra package:\n", color = :(green))

    (values, vectors) = eigen(entryMatrix)

    println("Eigenvalues: ", values)
    println("Eigenvectors: ")
    display(vectors)
    println()
end

function test_results_of_jacobi_method(entryMatrix::Array{Float64, 2})
    (eigenValues, eigenVectors) = jacobiMethod(entryMatrix)

    printstyled("JACOBI METHOD:\n", color = :(green))

    println("Computed values:")
    println("Eigenvalue: ", eigenValues)
    println("Eigenvector: ")
    display(eigenVectors)
    println()
    
    verifyEigenvalues(entryMatrix, eigenValues)
    verifyEigenvectors(entryMatrix, eigenVectors)
end

function test_results_of_power_method(entryMatrix::Array{Float64, 2})
    (eigenValues, eigenVectors) = powerMethodWithDeflation(entryMatrix)

    printstyled("POWER METHOD WITH DEFLATION:\n", color = :(green))

    println("Computed values:")
    println("Eigenvalue: ", eigenValues)
    println("Eigenvector: ")
    display(eigenVectors)
    println()

    verifyEigenvalues(entryMatrix, eigenValues)
    verifyEigenvectors(entryMatrix, eigenVectors)
end

function test_results_of_qr_method(entryMatrix::Array{Float64, 2})
    (eigenValues, eigenVectors) = qrMethod(entryMatrix)

    printstyled("QR METHOD\n", color = :(green))

    println("Computed values:")
    println("Eigenvalues: ", eigenValues)
    println("Eigenvectors: ")
    display(eigenVectors)
    println()

    verifyEigenvalues(entryMatrix, eigenValues)
    verifyEigenvectors(entryMatrix, eigenVectors)
end

function tests()
    entryMatrix::Array{Float64, 2} = Array(Symmetric(rand(4, 4)))

    eigenvalues_from_LinearAlgebra(entryMatrix)
    println()
    test_results_of_jacobi_method(entryMatrix)
    println()
    test_results_of_power_method(entryMatrix)
    println()
    test_results_of_qr_method(entryMatrix)
end

tests()
