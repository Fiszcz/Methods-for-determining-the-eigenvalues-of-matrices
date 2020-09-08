using LinearAlgebra

function jacobiMethod(A::Symmetric{Float64,Array{Float64,2}})
    (n, n) = size(A)
    A_copy = copy(A)
    previous_diagonal_length = typemax(Float64)
    current_diagonal_length_squared = 0
    eigenvalues = Array{Float64, 1}(undef, n)
    eigenvectors = Matrix{Float64}(I, n, n)
    while abs(previous_diagonal_length - current_diagonal_length_squared) > 0.001
        (s, t) = max_element_above_diagonal(A_copy)
        (sint, cost) = get_sin_cos(A_copy, s, t)
        R = Matrix{Float64}(I, n, n)
        R[s, s] = cost
        R[t, t] = cost
        R[s, t] = -sint
        R[t, s] = sint
        A_copy = Symmetric(transpose(R) * A_copy * R)
        eigenvectors = eigenvectors * R
        
        previous_diagonal_length = current_diagonal_length_squared
        current_diagonal_length_squared = 0
        for i = 1:n
            eigenvalues[i] = A_copy[i, i]
            current_diagonal_length_squared += A_copy[i, i]^2
        end
    end

    return (sort!(eigenvalues), eigenvectors)
end

function max_element_above_diagonal(A::Symmetric{Float64,Array{Float64,2}})
    x::Int64 = 1
    y::Int64 = 2
    (n, n) = size(A)
    for c = 2:n, r = 1:c-1
        if abs(A[r, c]) > abs(A[x, y])
            x = r
            y = c
        end
    end
    return (x, y)
end

function get_sin_cos(A::Symmetric{Float64,Array{Float64,2}}, s::Int64, t::Int64)
    d = sqrt((A[s, s] - A[t, t])^2 + 4*A[s, t]^2)
    sin2t = 2*A[s, t] / d
    cos2t = (A[s, s] - A[t, t]) / d
    dt = sqrt(2 * (1 + cos2t))
    sint = abs(sin2t) / dt
    cost = sign(A[s, t]) * abs((1 + cos2t) / dt)
    return (sint, cost)
end
