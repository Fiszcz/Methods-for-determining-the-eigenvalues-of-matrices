include("../power_method/powerMethod.jl")
include("../qr_method/qrMethod.jl")
include("./exampleMatrices.jl")

using Pkg
Pkg.add("BenchmarkTools")
using BenchmarkTools

printstyled("POWER METHOD:\n", color = :(green))
printstyled("Matrix 5x5:\n", color = :(blue))
display(@benchmark powerMethod($matrix5x5))
printstyled("Matrix 25x25:\n", color = :(blue))
display(@benchmark powerMethod($matrix25x25))
printstyled("Matrix 75x75:\n", color = :(blue))
display(@benchmark powerMethod($matrix75x75))
printstyled("Matrix 150x150:\n", color = :(blue))
display(@benchmark powerMethod($matrix150x150))

printstyled("QR METHOD:\n", color = :(green))
printstyled("Matrix 5x5:\n", color = :(blue))
display(@benchmark qrMethod($matrix5x5))
printstyled("Matrix 25x25:\n", color = :(blue))
display(@benchmark qrMethod($matrix25x25))
printstyled("Matrix 75x75:\n", color = :(blue))
display(@benchmark qrMethod($matrix75x75))
printstyled("Matrix 150x150:\n", color = :(blue))
display(@benchmark qrMethod($matrix150x150))
