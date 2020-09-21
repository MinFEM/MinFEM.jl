using MinFEM

mesh = import_mesh("../meshes/rounded.msh")

L = asmLaplacian(mesh)

f(x) = x[1]^2 + x[2]^2
s = evaluateMeshFunction(mesh, f, region=union(mesh.Boundaries[1003].Nodes,
                                               mesh.Boundaries[1004].Nodes))

s = asmBoundarySource(mesh, s, union(mesh.Boundaries[1003].Edges,
                                     mesh.Boundaries[1004].Edges))

pde = PDESystem(A=L, b=s, bc=zeros(mesh.nnodes), DI=union(mesh.Boundaries[1001].Nodes,
                                                          mesh.Boundaries[1002].Nodes))
solve(pde)

vtkfile = open_vtk_file(mesh, "boundary_source.vtu")
write_point_data(vtkfile, pde.state, "y")
write_point_data(vtkfile, s, "s")
save_vtk_file(vtkfile)

