using MinFEM
using WriteVTK

mesh = import_mesh("../meshes/poisson.msh")

L = asmElasticity(mesh,1.0,1.0)
Mb = asmBoundaryMassMatrix(mesh, mesh.Boundaries[1001].Edges, qdim=2)

f(x) = [0.1,0.0]
s = evaluateMeshFunction(mesh, f, region=mesh.Boundaries[1001].Nodes, qdim=2)

boundary = mesh.Boundaries[1003].Nodes

pde = PDESystem(L, Mb*s, zeros(2*mesh.nnodes), boundary, 2)
solve(pde)

vtkfile = write_vtk_mesh(mesh, "elasticity.vtu")
vtk_point_data(vtkfile, [reshape(pde.state, 2, mesh.nnodes); zeros(1,mesh.nnodes)], "Y")
vtk_point_data(vtkfile, s, "S")
vtk_save(vtkfile)

