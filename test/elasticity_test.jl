using MinFEM

function solve_elasticity(file_name::String)

  mesh = import_mesh(file_name)

  L = asmElasticity(mesh,1.0,1.0)
  Mb = asmBoundaryMassMatrix(mesh, Edges=mesh.Boundaries[1001].Edges, qdim=2)

  f(x) = [x[1]; x[2]]
  s = evaluateMeshFunction(mesh, f, region=mesh.Boundaries[1001].Nodes, qdim=2)

  boundary = mesh.Boundaries[1003].Nodes

  pde = PDESystem(A=L, b=Mb*s, bc=zeros(2*mesh.nnodes), DI=boundary, qdim=2)
  solve(pde)

  return true

end

@test solve_elasticity("test_v4.msh")
