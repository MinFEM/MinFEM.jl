using MinFEM

function solve_poisson(file_name::String)

  mesh = import_mesh(file_name)

  L = asmLaplacian(mesh)
  M = asmMassMatrix(mesh)

  n=3
  m=2
  f(x) = sin(n*x[1]*pi)*sin(m*x[2]*pi)
  s = evaluateMeshFunction(mesh, f)

  boundary = union(mesh.Boundaries[1001].Nodes,
                   mesh.Boundaries[1002].Nodes,
                   mesh.Boundaries[1003].Nodes,
                   mesh.Boundaries[1004].Nodes)

  pde = PDESystem(A=L, b=M*s, bc=zeros(mesh.nnodes), DI=boundary)

  solve(pde)

  v = ((n*pi)^2 + (m*pi)^2)*pde.state - s
  return sqrt(v'*M*v)

end

@test solve_poisson("test_v2.msh") == solve_poisson("test_v4.msh")
@test solve_poisson("test_v2.msh") < 4e-3
