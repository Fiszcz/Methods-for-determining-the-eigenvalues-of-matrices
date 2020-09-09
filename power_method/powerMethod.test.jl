include("../utils/verifyEigenvalues.jl")
include("./powerMethod.jl")
using LinearAlgebra
using BenchmarkTools
import Random

Random.seed!(100)
A = Array(Symmetric(rand(5, 5) * 100))

(eigenvalues, eigenvectors) = powerMethodWithDeflation(A)

display(eigvals(A))

display(eigenvectors)

println(eigenvalues)

display(@benchmark powerMethodWithDeflation(A))
