using LinearAlgebra

A = [2.0 1.0; 1.0 2.0]

A

function eigenvalue(A::Array{Float64,2}, v::Vector{Float64})::Array{Float64}
    return A * v
end

function powerIteration(A::Array{Float64,2}, size::Int64)
    v = ones(size) / sqrt(size)

    ev = eigenvalue(A, v)

    evNew::Array{Float64} = []
    vNew::Array{Float64} = []

    while true
        Av = A * v

        vNew = Av / norm(Av)

        evNew = eigenvalue(A, vNew)

        areAllWithPrecision::Bool = true
        for i=1:size
            if abs(ev[i] - evNew[i]) > 0.001
                areAllWithPrecision = false
                break
            end
        end
        if areAllWithPrecision
            break
        end

        v = vNew
        ev = evNew
    end

    return (evNew, vNew)
end

println(powerIteration(A, 2))
