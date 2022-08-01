using MinFEM
using LinearAlgebra

function semilinear(mesh::Mesh, L::AbstractMatrix, M::AbstractMatrix,
                    s::AbstractVector, BoundaryIndices::Set{Int64}=[], tol=1e-10)

    y = zeros(mesh.nnodes)

    pde = PDESystem(A=L, b=M*s, bc=zeros(mesh.nnodes), DI=BoundaryIndices)

    res = Inf
    while res > tol
        pde.A = L + assemble_cubicderivativematrix(mesh, y)
        pde.b = -L*y + M*s - assemble_cubicterm(mesh, y)
        refresh!(pde)
        solve!(pde)

        y += pde.state
        res = norm(pde.state)
        println(res)
    end
    
    return y
end

mesh = import_mesh("../meshes/Zshaped.msh")

L = assemble_laplacian(mesh)
M = assemble_massmatrix(mesh)

# y = 3*sin(x[1]*pi)*sin(x[2]*pi)
f(x) = 3*2*pi^2*sin(x[1]*pi)*sin(x[2]*pi) + (3*sin(x[1]*pi)*sin(x[2]*pi))^3
s = evaluate_mesh_function(mesh, f)

boundary = select_boundary(mesh)
boundaryNodes = extract_nodes(boundary)

y = semilinear(mesh, L, M, s, boundaryNodes);

write_to_vtk([y, s], mesh, ["y","s"], "semilinear")
