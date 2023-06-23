using MinFEM
using Test

tests = [
    "mesh_test",
    "quadrature_test",
    "evaluation_test",
    "properties_test",
    "deformation_test",
    "normal_test",
    "fem_test",
    "extension_test",
    "gradient_test",
    "poisson_test",
    "elasticity_test",
    "semilinear_test",
    "txtoutput_test",
    "vtkoutput_test"
]

for t in tests
    include("$(t).jl")
end
