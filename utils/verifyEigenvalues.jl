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

A = [2.0 3.0 5.0; 3.0 8.0 1.0; 5.0 1.0 3.0]
# First three element = correct, another three = wrong
eigenvalues = [10.628, 5.137, -2.766, 43, 112, -3]
# First element = correct, another one = wrong
eigenvectors = [[-0.760; 0.153; 0.632], [32.0; 64.0; 22.0]]

println("Verification for eigenvectors:")
verifyEigenvalues(A, eigenvalues)
println("\nVerification for eigenvectors:")
verifyEigenvectors(A, eigenvectors)
