include("../utils/verifyEigenvalues.jl")
include("./jacobiMethod.jl")
using LinearAlgebra
using BenchmarkTools
using InteractiveUtils
import Random

Random.seed!(100)
const A = Array(Symmetric(rand(5, 5) * 100))

(eigenvalues, eigenvectors) = jacobiMethod(A)

display(@code_warntype jacobiMethod(A))

display(eigvals(A))

display(eigenvectors)

println(eigenvalues)

Random.seed!(100)
const B = Array(Symmetric(rand(75, 75)))
display(@benchmark jacobiMethod(B))
