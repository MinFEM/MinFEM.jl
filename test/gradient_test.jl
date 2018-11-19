using MinFEM
using LinearAlgebra

function solve_gradient(file_name::String)

  mesh = import_mesh(file_name)

  G = asmGradient(mesh, qdim=2)

  f(x) = [x[1]+x[2]; x[1]-x[2]]
  s = evaluateMeshFunction(mesh, f, qdim=2)

  y = G*s;

  return y

end

@test norm(abs.(solve_gradient("test_v4.msh")).-1.0, Inf) < 1e-10
