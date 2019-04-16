using MinFEM
using WriteVTK

mesh = unit_square(100)

L = asmLaplacian(mesh)
M = asmMassMatrix(mesh)

n=3
m=2
f(x) = ((n*pi)^2 + (m*pi)^2) *sin(n*x[1]*pi)*sin(m*x[2]*pi)
s = evaluateMeshFunction(mesh, f)

boundary = union(mesh.Boundaries[1001].Nodes,
                  mesh.Boundaries[1002].Nodes,
                  mesh.Boundaries[1003].Nodes,
                  mesh.Boundaries[1004].Nodes)

pde = PDESystem(A=L, b=M*s, bc=zeros(mesh.nnodes), DI=boundary)
solve(pde)

vtkfile = write_vtk_mesh(mesh, "poisson.vtu")
vtk_point_data(vtkfile, pde.state, "Y")
vtk_point_data(vtkfile, s, "S")
vtk_save(vtkfile)

