using MinFEM
using LinearAlgebra

function test_pdesystem()
    mesh = unit_square(10)

    L = assemble_laplacian(mesh)
    M = assemble_massmatrix(mesh)

    f(x) = 0.0
    fh = evaluate_mesh_function(mesh, f)
    boundary = select_boundaries(mesh)

    pde = PDESystem(A=L, b=M*fh)
    assemble!(pde)

    pde.SystemMatrix != L && return false

    pde.bc = zeros(mesh.nnodes)
    pde.DI = extract_nodes(boundary)
    refresh!(pde)

    solve!(pde)
    (norm(pde.state) >= 1e-10) && return false
    
    return true
end

@test test_pdesystem()
