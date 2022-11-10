using MinFEM

function parabolic(mesh::Mesh, boundaryIndices::Set{Int64},
                    L::AbstractMatrix, M::AbstractMatrix,
                    f::AbstractVector, u::AbstractVector,
                    T::Float64, tsteps::Int; theta=1.0)

    # time increment assuming that t0=0
    dt = (T-0)/tsteps

    path = "parabolic_problem/"
    mkpath(path)
    write_to_vtk(u, mesh, "u", path*"parabolic_"*lpad(string(0), 3, '0')*".vtu")

    for i = 1:tsteps
        # we have to solve the following equation depending on theta
        # with the prescribed Dirichlet conditions
        # M*(u_new - u)/dt + L(theta*u_new + (1.0-theta)*u) = M*f
        pde = PDESystem(A=(M + theta*dt*L), b=(M - (1.0-theta)*dt*L)*u + dt*M*f,
                        bc=zeros(mesh.nnodes), DI=boundaryIndices)
        solve!(pde)

        u = copy(pde.state)

        write_to_vtk(u, mesh, "u", path*"parabolic_"*lpad(string(i), 3, '0')*".vtu")
    end

end

mesh = import_mesh("../meshes/Lshaped.msh")

boundary = select_boundaries(mesh)
boundary_nodes = extract_nodes(boundary)

source(x) = 1.0
f = evaluate_mesh_function(mesh, source)

initial_condition(x) = 0.0
u = evaluate_mesh_function(mesh, initial_condition)

L = assemble_laplacian(mesh)
M = assemble_massmatrix(mesh)

parabolic(mesh, boundary_nodes, L, M, f, u, 1.5, 100, theta=1.0)
