using MinFEM
using LinearAlgebra

function semilinear(mesh::Mesh, L::AbstractMatrix, M::AbstractMatrix,
                    s::AbstractVector, BoundaryIndices::Set{Int64}=[], tol=1e-10)

  y = zeros(mesh.nnodes)

  pde = PDESystem(A=L, b=M*s, bc=zeros(mesh.nnodes), DI=BoundaryIndices)

  res = Inf
  while res > tol
    pde.A = L + asmCubicDerivativeMatrix(mesh, y)
    pde.b = -L*y + M*s - asmCubicTerm(mesh, y)
    refresh(pde)
    solve(pde)

    y += pde.state
    res = norm(pde.state)
  end
  return y
end

function solve_semilinear(file_name::String)

  mesh = import_mesh(file_name)

  L = asmLaplacian(mesh)
  M = asmMassMatrix(mesh)

  f(x) = 100.0*x[1]*x[2]
  s = evaluateMeshFunction(mesh, f)

  boundary = union!(mesh.Boundaries[1001].Nodes,
                    mesh.Boundaries[1002].Nodes,
                    mesh.Boundaries[1003].Nodes,
                    mesh.Boundaries[1004].Nodes);

  y = semilinear(mesh, L, M, s, boundary);
  return true
end

@test solve_semilinear("test_v4.msh")
