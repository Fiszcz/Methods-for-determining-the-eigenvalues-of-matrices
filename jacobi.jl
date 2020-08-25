using LinearAlgebra

precisionOfAlgorithm = 0.0001

function Jacobi(A::Array{Float64, 2}, size::Int64)
    identityMatrix = Matrix{Float64}(I, 3, 3)
    s = copy(identityMatrix)

    flag = true
    while flag
        flag = false
    maxElementI = 1
    maxElementJ = 2
    maxValue = abs(A[1, 2])
    for j=1:size, i=1:size
        if i !== j && maxValue < abs(A[i, j])
            maxValue = abs(A[i, j])
            maxElementI = i
            maxElementJ = j
        end
    end

    θ = -π/4
    if A[maxElementI, maxElementI] === A[maxElementJ, maxElementJ] && A[maxElementI, maxElementJ] > 0
        θ = π/4
    else
        θ = 0.5 * atan(2 * A[maxElementI, maxElementJ] / (A[maxElementI, maxElementI] - A[maxElementJ, maxElementJ]))
    end

    s1 = copy(identityMatrix)
    s1t = copy(identityMatrix)

    cosθ = cos(θ)
    s1[maxElementI, maxElementI] = cosθ
    s1[maxElementJ, maxElementJ] = cosθ
    sinθ = sin(θ)
    s1[maxElementJ, maxElementI] = sinθ
    s1[maxElementI, maxElementJ] = -sinθ

    s1t = copy(s1)
    s1t[maxElementI, maxElementJ] = s1[maxElementJ, maxElementI]
    s1t[maxElementJ, maxElementI] = s1[maxElementI, maxElementJ]

    temp = zeros(size, size)
    for i=1:size, j=1:size
        sum = 0.0
        for p=1:size
            sum += s1t[i, p] * A[p, j]
        end
        temp[i, j] = sum
    end

    for i=1:size, j=1:size
        sum = 0.0
        for p=1:size
            sum += temp[i, p] * s1[p, j]
        end
        A[i, j] = sum
    end

    for i=1:size, j=1:size
        sum = 0.0
        for p=1:size
            sum += s[i, p] * s1[p, j]
        end
        temp[i, j] = sum
    end

    s = temp

    for i=1:size, j=1:size
        if i !== j && abs(A[i, j]) > precisionOfAlgorithm
            flag = true
        end
    end

end

eigenvalues = diag(A)
eigenvectors::Vector{Vector{Float64}} = []

for j=1:size
     push!(eigenvectors, s[:,j])
end

return (eigenvalues, eigenvectors)

end

A = [2.0 1.0; 1.0 2.0]
B = [5.0 3.0; 3.0 5.0]
C = [1.0 2.0 3.0; 2.0 4.0 2.0; 3.0 2.0 5.0]

println(Jacobi(C, 3))
