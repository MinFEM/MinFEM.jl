using MinFEM

function test_evaluation()
    id1(x) = x[1]
    id2(x) = [x[1], x[2]]
    id3(x) = [x[1], x[2], x[3]]

    mesh = import_mesh("test_line_v4.msh")
    dom = select_domains(mesh, 10001)
    for order = 1:9
        points = quadrature_points(mesh, order)
        xle = length(quadrature_points(mesh.d, order))

        vm = evaluate_mesh_function(mesh, id1)
        vq = evaluate_quadrature_function(mesh, id1, order=order)

        vq != evaluate_quadrature_function(mesh, id1, dom, order=order) && return false

        E = assemble_basismatrix(mesh, order = order)
        vmq = E * vm
        t = Array{Array{Float64,1},1}(undef,mesh.nelems*xle)
        for i = 1:mesh.nelems*xle
            t[i] = [vmq[i]]
        end

        for (i,x) in enumerate(points)
            if sum(abs.(t[i] .- x)) > 1e-14
                return false
            end
        end

        for (i,x) in enumerate(vq)
            if sum(abs.(vmq[i] .- x)) > 1e-14
                return false
            end
        end
    end
    
    mesh = import_mesh("test_square_v4.msh")
    dom = select_domains(mesh, 10001)
    for order = 1:8
        points = quadrature_points(mesh, order)
        xle = length(quadrature_points(mesh.d, order))

        vm = evaluate_mesh_function(mesh, id2, qdim=2)
        vq = evaluate_quadrature_function(mesh, id2, order=order, qdim=2)

        vq != evaluate_quadrature_function(mesh, id2, dom, order=order, qdim = 2) && 
            return false

        E = assemble_basismatrix(mesh, order = order, qdim=2)
        vmq = E * vm
        t = Array{Array{Float64,1},1}(undef,mesh.nelems*xle)
        for i = 1:mesh.nelems*xle
            t[i] = vmq[2*(i-1)+1 : 2*i]
        end

        for (i,x) in enumerate(points)
            if sum(abs.(t[i] .- x)) > 1e-14
                return false
            end
        end

        for (i,x) in enumerate(vq)
            if sum(abs.(vmq[i] .- x)) > 1e-14
                return false
            end
        end
    end

    mesh = import_mesh("test_cube_v4.msh")
    dom = select_domains(mesh, 10001)
    for order = 1:7
        points = quadrature_points(mesh, order)
        xle = length(quadrature_points(mesh.d, order))

        vm = evaluate_mesh_function(mesh, id3, qdim=3)
        vq = evaluate_quadrature_function(mesh, id3, order=order, qdim=3)

        vq != evaluate_quadrature_function(mesh, id3, dom, order=order, qdim = 3) && 
            return false

        E = assemble_basismatrix(mesh, order = order, qdim=3)
        vmq = E * vm
        t = Array{Array{Float64,1},1}(undef,mesh.nelems*xle)
        for i = 1:mesh.nelems*xle
            t[i] = vmq[3*(i-1)+1 : 3*i]
        end

        for (i,x) in enumerate(points)
            if sum(abs.(t[i] .- x)) > 1e-14
                return false
            end
        end

        for (i,x) in enumerate(vq)
            if sum(abs.(vmq[i] .- x)) > 1e-14
                return false
            end
        end
    end

    return true
end

function test_evaluation_boundary()
    id1(x) = x[1]
    id2(x) = [x[1], x[2]]
    id3(x) = [x[1], x[2], x[3]]

    mesh = import_mesh("test_line_v4.msh")
    boundary = select_boundaries(mesh)
    for order = 1:9
        points = quadrature_points_boundary(mesh, order)
        xle = length(quadrature_points_boundary(mesh.d, order))

        (
            points != quadrature_points_boundary(
                mesh,
                order,
                boundaryElements=extract_elements(boundary)
            )
        ) && return false

        vm = evaluate_mesh_function(mesh, id1, boundary)
        vq = evaluate_quadrature_function_boundary(mesh, id1, boundary, order=order)

        E = assemble_basismatrix_boundary(
            mesh,
            order = order,
            boundaryElements=extract_elements(boundary)
        )
        vmq = E * vm
        t = Array{Array{Float64,1},1}(undef,mesh.nboundelems*xle)
        for i = 1:mesh.nboundelems*xle
            t[i] = [vmq[i]]
        end

        for (i,x) in enumerate(points)
            if sum(abs.(t[i] .- x)) > 1e-14
                return false
            end
        end

        for (i,x) in enumerate(vq)
            if sum(abs.(vmq[i] .- x)) > 1e-14
                return false
            end
        end
    end
    boundary = select_boundaries(mesh, 1002)
    for order = 1:9
        xle = length(quadrature_points_boundary(mesh.d, order))

        vm = evaluate_mesh_function(mesh, id1, boundary)
        vq = evaluate_quadrature_function_boundary(mesh, id1, boundary, order=order)

        E = assemble_basismatrix_boundary(
            mesh,
            order = order,
            boundaryElements=extract_elements(boundary)
        )
        vmq = E * vm
        t = Array{Array{Float64,1},1}(undef,mesh.nboundelems*xle)
        for i = 1:mesh.nboundelems*xle
            t[i] = [vmq[i]]
        end

        for (i,x) in enumerate(vq)
            if sum(abs.(vmq[i] .- x)) > 1e-14
                return false
            end
        end
    end

    mesh = import_mesh("test_square_v4.msh")
    boundary = select_boundaries(mesh)
    for order = 1:8
        points = quadrature_points_boundary(mesh, order)
        xle = length(quadrature_points_boundary(mesh.d, order))

        (
            points != quadrature_points_boundary(
                mesh,
                order,
                boundaryElements=extract_elements(boundary)
            )
        ) && return false

        vm = evaluate_mesh_function(mesh, id2, boundary, qdim=2)
        vq = evaluate_quadrature_function_boundary(mesh, id2, boundary, order=order, qdim=2)

        E = assemble_basismatrix_boundary(
            mesh,
            order = order,
            boundaryElements=extract_elements(boundary),
            qdim=2
        )
        vmq = E * vm
        t = Array{Array{Float64,1},1}(undef,mesh.nboundelems*xle)
        for i = 1:mesh.nboundelems*xle
            t[i] = vmq[2*(i-1)+1 : 2*i]
        end

        for (i,x) in enumerate(points)
            if sum(abs.(t[i] .- x)) > 1e-14
                return false
            end
        end

        for (i,x) in enumerate(vq)
            if sum(abs.(vmq[i] .- x)) > 1e-14
                return false
            end
        end
    end
    boundary = select_boundaries(mesh, 1001, 1003)
    for order = 1:8
        xle = length(quadrature_points_boundary(mesh.d, order))

        vm = evaluate_mesh_function(mesh, id2, boundary, qdim=2)
        vq = evaluate_quadrature_function_boundary(mesh, id2, boundary, order=order, qdim=2)

        E = assemble_basismatrix_boundary(
            mesh,
            order = order,
            boundaryElements=extract_elements(boundary), 
            qdim=2
        )
        vmq = E * vm
        t = Array{Array{Float64,1},1}(undef,mesh.nboundelems*xle)
        for i = 1:mesh.nboundelems*xle
            t[i] = vmq[2*(i-1)+1 : 2*i]
        end

        for (i,x) in enumerate(vq)
            if sum(abs.(vmq[i] .- x)) > 1e-14
                return false
            end
        end
    end

    mesh = import_mesh("test_cube_v4.msh")
    boundary = select_boundaries(mesh)
    for order = 1:7
        points = quadrature_points_boundary(mesh, order)
        xle = length(quadrature_points_boundary(mesh.d, order))

        (
            points != quadrature_points_boundary(
                mesh,
                order,
                boundaryElements=extract_elements(boundary)
            )
        ) && return false

        vm = evaluate_mesh_function(mesh, id3, boundary, qdim=3)
        vq = evaluate_quadrature_function_boundary(mesh, id3, boundary, order=order, qdim=3)

        E = assemble_basismatrix_boundary(
            mesh,
            order = order,
            boundaryElements=extract_elements(boundary),
            qdim=3
        )
        vmq = E * vm
        t = Array{Array{Float64,1},1}(undef,mesh.nboundelems*xle)
        for i = 1:mesh.nboundelems*xle
            t[i] = vmq[3*(i-1)+1 : 3*i]
        end

        for (i,x) in enumerate(points)
            if sum(abs.(t[i] .- x)) > 1e-14
                return false
            end
        end

        for (i,x) in enumerate(vq)
            if sum(abs.(vmq[i] .- x)) > 1e-14
                return false
            end
        end
    end
    boundary = select_boundaries(mesh, 1001, 1003, 1004)
    for order = 1:7
        xle = length(quadrature_points_boundary(mesh.d, order))

        vm = evaluate_mesh_function(mesh, id3, boundary, qdim=3)
        vq = evaluate_quadrature_function_boundary(mesh, id3, boundary, order=order, qdim=3)

        E = assemble_basismatrix_boundary(
            mesh,
            order = order,
            boundaryElements=extract_elements(boundary), 
            qdim=3
        )
        vmq = E * vm
        t = Array{Array{Float64,1},1}(undef,mesh.nboundelems*xle)
        for i = 1:mesh.nboundelems*xle
            t[i] = vmq[3*(i-1)+1 : 3*i]
        end

        for (i,x) in enumerate(vq)
            if sum(abs.(vmq[i] .- x)) > 1e-14
                return false
            end
        end
    end
    
    return true
end

@test test_evaluation()
@test test_evaluation_boundary()
