A = [2.0 3.0 5.0; 3.0 8.0 1.0; 5.0 1.0 3.0]
# First three element = correct, another three = wrong
eigenvalues = [10.628, 5.137, -2.766, 43, 112, -3]
# First element = correct, another one = wrong
eigenvectors = [[-0.760; 0.153; 0.632], [32.0; 64.0; 22.0]]

println("Verification for eigenvectors:")
verifyEigenvalues(A, eigenvalues)
println("\nVerification for eigenvectors:")
verifyEigenvectors(A, eigenvectors)
