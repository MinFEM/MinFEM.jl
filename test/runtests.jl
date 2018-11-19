using MinFEM
using Test

tests = ["poisson_test", "elasticity_test", "semilinear_test", "gradient_test"]
for t in tests
  include("$(t).jl")
end
