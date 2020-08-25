using LinearAlgebra

ε = 2e-6

function householder(matrix::Array{Float64, 2}, size::Int64)
    for k=1:(size-2)
        if abs(matrix[k, k+1]) < ε
            continue
        end

        column = [zeros(k)..., matrix[k+1:size, k]...]

        σ = sign(column[k+1]) * hypot(column...)
        column[k+1] += σ
        norm = sqrt(2σ * column[k + 1])

        map!(x -> x / norm, column, column)

        sdx = zeros(size)
        sdy = zeros(size)
        for i=1:size, j=1:size
            sdx[i] += column[j] * matrix[j, i]
            sdy[i] += column[j] * matrix[i, j]
        end

        hdx = 0.0
        hdy = 0.0
        for i=1:size
            hdx += sdx[i] * column[i]
            hdy += sdy[i] * column[i]
        end

        dx = Vector{Float64}(undef, size)
        dy = Vector{Float64}(undef, size)
        for i=1:size
            dx[i] = 2(sdx[i] - hdx * column[i])
            dy[i] = 2(sdy[i] - hdy * column[i])
        end

        for i=1:size, j=1:size
            matrix[j, i] -= column[i] * dy[j] + column[j] * dx[i]
        end
    end

    return matrix
end

function QR(matrix::Array{Float64, 2}, size::Int64)
    identityMatrix = Matrix{Float64}(I, size, size)
    k = size

    while k >= 2
        if (abs(matrix[k - 1, k]) < ε)
            k -= 1
            continue
        end

        a = matrix[k, k]
        b = matrix[k, k - 1]
        c = matrix[k - 1, k]
        d = matrix[k - 1, k - 1]

        tr = a * d
        det = a * d - b * c
        disc = sqrt(tr^2 - 4det)
        l1 = (tr + disc) / 2
        l2 = (tr - disc) / 2
        mu = a - (abs(l1) < abs(l2) ? l1 : l2)

        matrix -= mu * I

        q = Matrix{Float64}(I, size, size)
        for i=1:k-1
            a1 = matrix[i, i]
            a2 = matrix[i + 1, i]
            base = hypot(a1, a2)
            cos = base < ε ? 0 : a1 / base
            sin = base < ε ? 0 : a2 / base

            matrix[i, i] = base
            matrix[i + 1, i] = 0
            for j=i+1:k
                e1 = matrix[i, j]
                e2 = matrix[i + 1, j]
                matrix[i, j] = e1 * cos + e2 * sin
                matrix[i + 1, j] = e2 * cos - e1 * sin
            end

            for j=1:k
                e1 = q[j, i]
                e2 = q[j, i + 1]
                q[j, i] = e1 * cos + e2 * sin
                q[j, i + 1] = e2 * cos - e1 * sin
            end
        end
        println(q)

        matrix = matrix * q
        matrix += mu * I
    end

    return matrix
end

function QRForEigenvalues(matrix::Array{Float64, 2})
    (n, n) = size(matrix)
    println(householder(matrix, n))
    println(diag(QR(matrix, n)))
end

A = [2.0 1.0; 1.0 2.0]
B = [5.0 3.0; 3.0 5.0]
C = [1.0 2.0 3.0; 2.0 4.0 2.0; 3.0 2.0 5.0]

QRForEigenvalues(C)
