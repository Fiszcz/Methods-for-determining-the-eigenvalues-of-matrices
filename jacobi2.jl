using LinearAlgebra

function jacobi(A::Array{Float64, 2})
    (n, n) = size(A)
    result = copy(A)
    for i = 1:70 # TODO: better convergence checking method needed (probably )
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
    end

    eigenvalues = Array{Float64, 1}(undef, n)
    for i = 1:n
        eigenvalues[i] = result[i, i]
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

function main()
    A::Array{Float64, 2} = [
        12 6 -6 1
        6 16 2 100
        -6 2 16 1
        1 100 1 1]

    println("A: $A")
    println("Correct eigenvalues: $(eigvals(A))")
    println("Jacobi eigenvalues: $(jacobi(A))")
end

main()
