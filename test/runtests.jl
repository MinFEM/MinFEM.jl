using MinFEMDev
using Test

tests = ["quadrature_test", "mesh_test", "properties_test", "deformation_test", "normal_test", 
            "norm_test", "gradient_test", "poisson_test", "elasticity_test", "semilinear_test",
            "txtoutput_test", "vtkoutput_test"]

for t in tests
    include("$(t).jl")
end
