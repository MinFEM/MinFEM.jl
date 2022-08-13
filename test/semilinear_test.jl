using MinFEM
using LinearAlgebra

function test_semilinear(file_name::String)

    mesh = import_mesh(file_name)

    L = assemble_laplacian(mesh)
    M = assemble_massmatrix(mesh)

    f(x) = 100.0 * x[1] * x[2]
    s = evaluate_mesh_function(mesh, f)

    boundary = select_boundaries(mesh)

    pde = PDESystem(A=L, b=M*s, bc=zeros(mesh.nnodes), DI=extract_nodes(boundary))

    y = zeros(mesh.nnodes)
    res = Inf
    for it = 1:6
        pde.A = L + assemble_cubicderivativematrix(mesh, y)
        pde.b = -L*y + M*s - assemble_cubicterm(mesh, y)
        refresh!(pde)
        solve!(pde)

        y += pde.state
        res = norm(pde.state)
        
        if res <= 1e-10
            return true
        end
    end

    return false
end

@test test_semilinear("test_square_v4.msh")
