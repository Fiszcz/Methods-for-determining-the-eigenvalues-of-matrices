function verifyEigenvalues(entryMatrix::Matrix{Float64}, eigenvalues::Vector{Float64})
    for eigenvalue in eigenvalues
        determinant = det(entryMatrix - eigenvalue*I)
        if determinant < 0.05 && determinant > -0.05
            println("Correct eigenvalue: ", eigenvalue)
        else
            println("Wrong eigenvalue: ", eigenvalue)
            println("Determinant: ", determinant)
        end
    end
end

function verifyEigenvectors(entryMatrix::Matrix{Float64}, eigenvectors::Array{Vector{Float64}})
    eigenvalues = eigvals(A)
    for eigenvector in eigenvectors
        if any(eigenvalue -> all(element -> element > -0.05 && element < 0.05, (entryMatrix - eigenvalue*I)*eigenvector), eigenvalues)
            println("Correct eigenvector: ", eigenvector)
        else
            println("Wrong eigenvector: ", eigenvector)
        end
    end
end
