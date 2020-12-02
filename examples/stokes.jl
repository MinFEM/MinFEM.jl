using MinFEM

mesh = import_mesh("../meshes/stokes.msh");
L = asmVectorLaplacian(mesh, qdim=2);
B = asmIncompressibility(mesh);
C = MinFEM.asmStokesStabilization(mesh);
alpha = 0.05;

A = [L B; B' -alpha*C];

boundingBox = MinFEM.getMeshBoundingBox(mesh)
f(x) = [0.5*cos(2*pi*x[2]/(boundingBox[1][2]-boundingBox[2][2]))+0.5,0.0]
bc = evaluateMeshFunction(mesh, f, region=mesh.Boundaries[1003].Nodes,
                                   qdim=2)
                                   
rhs = zeros(3*mesh.nnodes)
asmDirichletCondition(A, union(mesh.Boundaries[1001].Nodes,
                               mesh.Boundaries[1002].Nodes,
                               mesh.Boundaries[1003].Nodes),
                      rhs=rhs, bc=bc, qdim=2)

x = A\rhs;

vtkfile = open_vtk_file(mesh, "stokes.vtu")
write_point_data(vtkfile, [reshape(x[1:2*mesh.nnodes], 2, mesh.nnodes); zeros(1,mesh.nnodes)], "v")
write_point_data(vtkfile, x[2*mesh.nnodes+1:end], "p")
save_vtk_file(vtkfile)
