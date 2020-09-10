using LinearAlgebra

function jacobiMethod(A::Array{Float64,2})
    (n, n) = size(A)
    A_copy = copy(A)

    eigenvalues = Array{Float64, 1}(undef, n)
    eigenvectors = Array{Float64, 2}(I, n, n)

    previous_diagonal_length::Float64 = 10.0
    current_diagonal_length_squared::Float64 = eps(Float64)
    while abs((previous_diagonal_length - current_diagonal_length_squared)/previous_diagonal_length) > 0.001
        (s, t) = max_element_above_diagonal(A_copy)
        (sint, cost) = get_sin_cos(A_copy, s, t)
        R = Matrix{Float64}(I, n, n)
        R[s, s] = cost
        R[t, t] = cost
        R[s, t] = -sint
        R[t, s] = sint
        A_copy = R' * A_copy * R
        eigenvectors *= R

        previous_diagonal_length = current_diagonal_length_squared
        current_diagonal_length_squared = sum(diag(A_copy).^2)
    end

    return (diag(A_copy), eigenvectors)
end

function max_element_above_diagonal(A::Array{Float64,2})
    row_with_max::Int64 = 1
    column_with_max::Int64 = 2
    max_value = abs(A[1, 2])
    (n, n) = size(A)
    for column = 3:n, row = 1:column-1
        if abs(A[row, column]) > max_value
            row_with_max = row
            column_with_max = column
            max_value = abs(A[row, column])
        end
    end
    return (row_with_max, column_with_max)
end

function get_sin_cos(A::Array{Float64,2}, s::Int64, t::Int64)
    d = sqrt((A[s, s] - A[t, t])^2 + 4*A[s, t]^2)
    sin2t = 2*A[s, t] / d
    cos2t = (A[s, s] - A[t, t]) / d
    dt = sqrt(2 * (1 + cos2t))
    sint = abs(sin2t) / dt
    cost = sign(A[s, t]) * abs((1 + cos2t) / dt)
    return (sint, cost)
end
