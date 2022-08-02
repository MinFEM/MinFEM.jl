using MinFEMDev

function test_elasticity(fileName::String)
    mesh = import_mesh(fileName)

    forcingBoundary = select_boundaries(mesh, 1001)
    fixedBoundary = select_boundaries(mesh, 1003)

    L = assemble_elasticity(mesh, 1.0, 1.0)
    Mb = assemble_massmatrix_boundary(mesh, boundaryElements=extract_elements(forcingBoundary), qdim=2)

    f(x) = [x[1]; x[2]]
    s = evaluate_mesh_function(mesh, f, region=extract_nodes(forcingBoundary), qdim=2)

    pde = PDESystem(A=L, b=Mb*s, bc=zeros(2*mesh.nnodes), DI=extract_nodes(fixedBoundary), qdim=2)
    solve!(pde)

    return true
end

@test test_elasticity("test_square_v4.msh")