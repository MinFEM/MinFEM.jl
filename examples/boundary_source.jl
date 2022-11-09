using MinFEM

mesh = import_mesh("../meshes/rounded.msh")

freeBoundary = select_boundaries(mesh, 1003, 1004)
freeElements = extract_elements(freeBoundary)
dirichletBoundary = select_boundaries(mesh, 1001, 1002)
dirichletNodes = extract_nodes(dirichletBoundary)

L = assemble_laplacian(mesh)
M = assemble_massmatrix_boundary(mesh, boundaryElements=freeElements)

f(x) = x[1]^2 + x[2]^2
s = evaluate_mesh_function(mesh, f, freeBoundary)

pde = PDESystem(A=L, b=M*s, bc=zeros(mesh.nnodes), DI=dirichletNodes)
solve!(pde)

write_to_vtk([pde.state, s], mesh, ["Y","S"], "boundary_source")
