include("../utils/verifyEigenvalues.jl")
include("./powerMethod.jl")
using LinearAlgebra
using BenchmarkTools
using InteractiveUtils
import Random

Random.seed!(100)
const A = Array(Symmetric(rand(5, 5) * 100))

(eigenvalues, eigenvectors) = powerMethodWithDeflation(A)

display(@code_warntype powerMethodWithDeflation(A))

display(eigvals(A))

display(eigenvectors)

println(eigenvalues)

Random.seed!(100)
const B = Array(Symmetric(rand(30, 30)))
display(@benchmark powerMethodWithDeflation(B))
