using LinearAlgebra

function jacobi2(A::Array{Float64, 2})
    (n, n) = size(A)
    result = copy(A)
    previous_diagonal_length = typemax(Float64)
    current_diagonal_length_squared = 0
    eigenvalues = Array{Float64, 1}(undef, n)
    while abs(previous_diagonal_length - current_diagonal_length_squared) > 0.001
        (s, t) = max_element_above_diagonal(result)
        d = sqrt((result[s, s] - result[t, t])^2 + 4*result[s, t]^2)
        sin2t = 2*result[s, t] / d
        cos2t = (result[s, s] - result[t, t]) / d
        dt = sqrt(2 * (1 + cos2t))
        sint = abs(sin2t) / dt
        cost = sign(result[s, t]) * abs((1 + cos2t) / dt)

        R = Matrix{Float64}(I, n, n)
        R[s, s] = cost
        R[t, t] = cost
        R[s, t] = -sint
        R[t, s] = sint
        result = transpose(R) * result * R
        
        previous_diagonal_length = current_diagonal_length_squared
        current_diagonal_length_squared = 0
        for i = 1:n
            eigenvalues[i] = result[i, i]
            current_diagonal_length_squared += result[i, i]^2
        end
    end

    return sort!(eigenvalues)
end

function max_element_above_diagonal(A::Array{Float64, 2})
    x = 1
    y = 2
    (n, n) = size(A)
    for c = 2:n, r = 1:c-1
        if abs(A[r, c]) > abs(A[x, y])
            x = r
            y = c
        end
    end
    return (x, y)
end

function jacobi(A::Array{Float64, 2})
    (n, n) = size(A)
    identityMatrix = Matrix{Float64}(I, n, n)
    s = copy(identityMatrix)
    precisionOfAlgorithm = 0.0001

    flag = true
    while flag
        flag = false

        maxElementI = 1
        maxElementJ = 2
        maxValue = abs(A[1, 2])
        for j=1:n, i=1:n
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

        temp = zeros(n, n)
        for i=1:n, j=1:n
            sum = 0.0
            for p=1:n
                sum += s1t[i, p] * A[p, j]
            end
            temp[i, j] = sum
        end

        for i=1:n, j=1:n
            sum = 0.0
            for p=1:n
                sum += temp[i, p] * s1[p, j]
            end
            A[i, j] = sum
        end

        for i=1:n, j=1:n
            sum = 0.0
            for p=1:n
                sum += s[i, p] * s1[p, j]
            end
            temp[i, j] = sum
        end

        s = temp

        for i=1:n, j=1:n
            if i !== j && abs(A[i, j]) > precisionOfAlgorithm
                flag = true
            end
        end

    end

    eigenvalues = diag(A)
    eigenvectors::Vector{Vector{Float64}} = []

    for j=1:n
         push!(eigenvectors, s[:,j])
    end

    return (sort!(eigenvalues), eigenvectors)

end

function main()
    A::Array{Float64, 2} = [
        12 6 -6 1
        6 16 2 100
        -6 2 16 1
        1 100 1 1]

    println("A: $A")
    println("Correct eigenvalues: $(@time (eigvals(A), eigvecs(A)))")
    println("Jacobi2 eigenvalues: $(@time jacobi2(A))")
    println("Jacobi eigenvalues: $(@time jacobi(A))")
end

main()
