using MinFEM

mesh = import_mesh("../meshes/square.msh")

forcingBoundary = select_boundaries(mesh, 1001)
fixedBoundary = select_boundaries(mesh, 1003)

L = assemble_elasticity(mesh, 1.0, 1.0)
Mb = assemble_massmatrix_boundary(mesh, boundaryElements=extract_elements(forcingBoundary), qdim=2)

f(x) = [0.1, 0.0]
s = evaluate_mesh_function(mesh, f, region=extract_nodes(forcingBoundary), qdim=2)

pde = PDESystem(A=L, b=Mb*s, bc=zeros(2*mesh.nnodes), DI=extract_nodes(fixedBoundary), qdim=2)
solve!(pde)

write_to_vtk([pde.state, s], mesh, ["Y","S"], "elasticity", qdim=2)
