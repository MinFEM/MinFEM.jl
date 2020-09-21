using MinFEM
#using WriteVTK

mesh = import_mesh("../meshes/poisson.msh")

L = asmElasticity(mesh,1.0,1.0)
Mb = asmBoundaryMassMatrix(mesh, Edges=mesh.Boundaries[1001].Edges, qdim=2)

f(x) = [0.1,0.0]
s = evaluateMeshFunction(mesh, f, region=mesh.Boundaries[1001].Nodes, qdim=2)

boundary = mesh.Boundaries[1003].Nodes

pde = PDESystem(A=L, b=Mb*s, bc=zeros(2*mesh.nnodes), DI=boundary, qdim=2)
solve(pde)

vtkfile = open_vtk_file(mesh, "elasticity.vtu")
write_point_data(vtkfile, [reshape(pde.state, 2, mesh.nnodes); zeros(1,mesh.nnodes)], "Y")
write_point_data(vtkfile, s, "S")
save_vtk_file(vtkfile)

