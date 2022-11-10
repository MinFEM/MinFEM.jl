using MinFEM
using LinearAlgebra

function semilinear(mesh::Mesh, L::AbstractMatrix, M::AbstractMatrix,
                    s::AbstractVector, boundaryIndices::Set{Int64};
                    tol::Float64=1e-10, maxIterations::Int64=10)

    y = zeros(mesh.nnodes)

    pde = PDESystem(A=L, b=M*s, bc=zeros(mesh.nnodes), DI=boundaryIndices)

    for i = 1:maxIterations
        pde.A = L + assemble_cubicderivativematrix(mesh, y)
        pde.b = -L*y + M*s - assemble_cubicterm(mesh, y)
        refresh!(pde)
        solve!(pde)

        y += pde.state
        res = norm(pde.state)

        if res < tol
            println("It. $i: $res < $tol")
            println("Semi-linear routine converged.")
            break
        else
            println("It. $i: $res ≥ $tol")

            if i == maxIterations
                println("Semi-linear routine failed.")
                println("Maximum number of iterations reached.")
            end
        end
    end
    
    return y
end

mesh = import_mesh("../meshes/Zshaped.msh")

L = assemble_laplacian(mesh)
M = assemble_massmatrix(mesh)

f(x) = 3*2*pi^2*sin(x[1]*pi)*sin(x[2]*pi) + (3*sin(x[1]*pi)*sin(x[2]*pi))^3
s = evaluate_mesh_function(mesh, f)

boundary = select_boundaries(mesh)
boundaryNodes = extract_nodes(boundary)

y = semilinear(mesh, L, M, s, boundaryNodes)

write_to_vtk([y, s], mesh, ["y","s"], "semilinear")
