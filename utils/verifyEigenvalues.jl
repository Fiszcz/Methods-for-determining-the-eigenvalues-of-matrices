function verifyEigenvalues(entryMatrix::Array{Float64, 2}, eigenvalues::Vector{Float64})
    for eigenvalue in eigenvalues
        determinant = det(entryMatrix - eigenvalue*I)
        if determinant < 0.05 && determinant > -0.05
            println("Correct eigenvalue: ", eigenvalue)
        else
            println("Wrong eigenvalue: ", eigenvalue)
        end
    end
end

function verifyEigenvectors(entryMatrix::Array{Float64, 2}, eigenvectors::Array{Float64, 2})
    eigenvalues = eigvals(entryMatrix)
    (n, n) = size(entryMatrix)
    for i = 1:n
        eigenvector = eigenvectors[:,i]
        if any(eigenvalue -> all(element -> element > -0.05 && element < 0.05, (entryMatrix - eigenvalue*I)*eigenvector), eigenvalues)
            println("Correct eigenvector: ", eigenvector)
        else
            println("Wrong eigenvector: ", eigenvector)
        end
    end
end
