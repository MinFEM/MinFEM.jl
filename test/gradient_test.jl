using MinFEMDev
using LinearAlgebra

function test_gradient()
    f1(x) = -x[1]
    g1(x) = -1
    f3(x) = [x[1]+x[2]+x[3]; x[1]-x[2]+x[3]; x[1]+x[2]-x[3]]


    mesh::Mesh = import_mesh("test_line_v4.msh")
    G = assemble_derivativematrix(mesh)
    fh = evaluate_mesh_function(mesh, f1)

    g = -ones(Float64, mesh.nelems)
    any(abs.(G*fh .- g) .> 1e-15) && return false


    mesh = import_mesh("test_square_v4.msh")
    G = assemble_derivativematrix(mesh)
    fh = evaluate_mesh_function(mesh, f1)

    g = zeros(Float64, mesh.d*mesh.nelems)
    for i = 1:mesh.nelems
        g[2*(i-1)+1:2*i] = [-1, 0]
    end

    any(abs.(G*fh .- g) .> 1e-14) && return false

    mesh = import_mesh("test_cube_v4.msh")
    G = assemble_derivativematrix(mesh)
    fh = evaluate_mesh_function(mesh, f1)

    g = zeros(Float64, mesh.d*mesh.nelems)
    for i = 1:mesh.nelems
        g[3*(i-1)+1:3*i] = [-1, 0, 0]
    end

    any(abs.(G*fh .- g) .> 1e-14) && return false

    G = assemble_derivativematrix(mesh, qdim=3)
    fh = evaluate_mesh_function(mesh, f3, qdim=3)

    g = zeros(Float64, 3*mesh.d*mesh.nelems)
    for i = 1:mesh.nelems
        g[3*3*(i-1)+1:3*3*i] = [1, 1, 1, 1, -1, 1, 1, 1, -1]
    end

    any(abs.(G*fh .- g) .> 1e-14) && return false

    return true
end

@test test_gradient()
