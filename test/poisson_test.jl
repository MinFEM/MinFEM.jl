using MinFEMDev

function test_poisson()
    n = 3
    m = 2
    f1(x) = sin(n*x[1]*pi)
    f2(x) = sin(n*x[1]*pi) * sin(m*x[2]*pi)

    mesh::Mesh = import_mesh("test_line_v4.msh")
    boundaryNodes = extract_nodes(select_boundaries(mesh))

    L = assemble_laplacian(mesh)
    M = assemble_massmatrix(mesh)

    s = evaluate_mesh_function(mesh, f1)
    rhs = M * s

    pde = PDESystem(A=L, b=rhs, bc=zeros(mesh.nnodes), DI=boundaryNodes)
    solve!(pde)
    
    assemble_dirichletcondition!(L, boundaryNodes, rhs=rhs, bc=zeros(mesh.nnodes))
    y = L \ rhs

    any(abs.(pde.state .- y) .> 1e-15) && return false

    v1 = (n*pi)^2 * pde.state - s
    sqrt(v1' * M * v1) > 6e-2 && return false

    v2 = (n*pi)^2 * y - s
    sqrt(v2' * M * v2) > 6e-2 && return false


    mesh = import_mesh("test_square_v4.msh")
    boundaryNodes = extract_nodes(select_boundaries(mesh))

    L = assemble_laplacian(mesh)
    M = assemble_massmatrix(mesh)

    s = evaluate_mesh_function(mesh, f2)
    rhs = M * s

    pde = PDESystem(A=L, b=rhs, bc=zeros(mesh.nnodes), DI=boundaryNodes)
    solve!(pde)
    
    assemble_dirichletcondition!(L, boundaryNodes, rhs=rhs, bc=zeros(mesh.nnodes))
    y = L \ rhs

    any(abs.(pde.state .- y) .> 1e-15) && return false

    v1 = ((n*pi)^2 + (m*pi)^2) * pde.state - s
    sqrt(v1' * M * v1) > 6e-2 && return false

    v2 = ((n*pi)^2 + (m*pi)^2) * y - s
    sqrt(v2' * M * v2) > 6e-2 && return false
    
    return true
end

@test test_poisson()
