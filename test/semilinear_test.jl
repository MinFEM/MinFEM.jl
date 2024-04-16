using MinFEMDev
using LinearAlgebra
using SparseArrays

function test_semilinear_newton()

    mesh = import_mesh("test_square_v4.msh")

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

function test_semilinear_lagrangenewton()

    L2norm(M,v) = sqrt(v' * M * v) 

    mesh = import_mesh("test_square_v4.msh")
  
    L = assemble_laplacian(mesh)
    M = assemble_massmatrix(mesh)
  
    alpha = 1e-6
  
    f(x) = 5.0*max(sin(3*x[1]*pi)*sin(2*x[2]*pi), 0)
    z = evaluate_mesh_function(mesh, f)
  
    n = mesh.nnodes
    Z = spzeros(n, n)
    u = zeros(n)
    y = zeros(n)
    p = zeros(n)
  
    for it = 1:10
        eyu = L*y + assemble_cubicterm(mesh, y) + M*y - M*u
        LD1 = L + assemble_cubicderivativematrix(mesh, y) + M
        D2 = assemble_cubicsecondderivativematrix(mesh, y, p)
  
        A = [M+D2          Z   -LD1;
                Z    alpha*M      M;
                -LD1          M      Z]
  
        rhs = [LD1*p - M*(y-z);  -M*p - alpha*M*u; eyu];
        delta = A\rhs;
  
        y += delta[0*n+1 : 1*n]
        u += delta[1*n+1 : 2*n]
        p += delta[2*n+1 : 3*n]
  
        if L2norm(M, delta[1*n+1:2*n]) < 1e-10
            return true
        end
    end
  
    return false
end

@test test_semilinear_newton()
@test test_semilinear_lagrangenewton()
