using MinFEM

mesh = unit_square(100)
L = assemble_laplacian(mesh)
M = assemble_massmatrix(mesh)

n=3
m=2
f(x) = ((n*pi)^2 + (m*pi)^2) * sin(n*x[1]*pi) * sin(m*x[2]*pi)
s = evaluate_mesh_function(mesh, f)

boundary = select_boundaries(mesh, 1001, 1002, 1003, 1004)
boundaryNodes = extract_nodes(boundary)

pde = PDESystem(A=L, b=M*s, bc=zeros(mesh.nnodes), DI=boundaryNodes)
solve!(pde)
write_to_vtk([pde.state, s], mesh, ["Y","S"], "poisson")
