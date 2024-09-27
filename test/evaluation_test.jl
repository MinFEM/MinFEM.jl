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

        for (i,x) in enumerate(mesh.Nodes)
            if sum(abs.(evaluate_function(vm, i) .- x)) > 1e-14 
                return false
            end
        end

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

        for (i,x) in enumerate(mesh.Nodes)
            if sum(abs.(evaluate_function(vm, i, qdim=2) .- x)) > 1e-14 
                return false
            end
        end

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

        for (i,x) in enumerate(mesh.Nodes)
            if sum(abs.(evaluate_function(vm, i, qdim=3) .- x)) > 1e-14 
                return false
            end
        end

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
    belems = extract_elements(boundary)

    vf = evaluate_mesh_function(mesh, id1)
    vm = evaluate_mesh_function(mesh, id1, boundary)

    if any(abs.(evaluate_function(vf, extract_nodes(boundary)) .- vm) .> 1e-14)
        return false
    end

    for order = 1:9
        points = quadrature_points_boundary(mesh, order)
        xle = length(quadrature_points_boundary(mesh.d, order))

        (
            points != quadrature_points_boundary(
                mesh,
                order,
                boundaryElements = belems
            )
        ) && return false

        vq = evaluate_quadrature_function_boundary(mesh, id1, boundary, order=order)

        E = assemble_basismatrix_boundary(
            mesh,
            order = order,
            boundaryElements = belems
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
    belems = extract_elements(boundary)

    vm = evaluate_mesh_function(mesh, id1, boundary)

    if any(abs.(evaluate_function(vf, extract_nodes(boundary)) .- vm) .> 1e-14)
        return false
    end

    for order = 1:9
        xle = length(quadrature_points_boundary(mesh.d, order))

        vq = evaluate_quadrature_function_boundary(mesh, id1, boundary, order=order)

        E = assemble_basismatrix_boundary(
            mesh,
            order = order,
            boundaryElements = belems
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
    belems = extract_elements(boundary)

    vf = evaluate_mesh_function(mesh, id2, qdim=2)
    vm = evaluate_mesh_function(mesh, id2, boundary, qdim=2)

    if any(abs.(evaluate_function(vf, extract_nodes(boundary), qdim=2) .- vm) .> 1e-14)
        return false
    end

    for order = 1:8
        points = quadrature_points_boundary(mesh, order)
        xle = length(quadrature_points_boundary(mesh.d, order))

        (
            points != quadrature_points_boundary(
                mesh,
                order,
                boundaryElements = belems
            )
        ) && return false

        vq = evaluate_quadrature_function_boundary(mesh, id2, boundary, order=order, qdim=2)

        E = assemble_basismatrix_boundary(
            mesh,
            order = order,
            boundaryElements = belems,
            qdim = 2
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
    belems = extract_elements(boundary)

    vm = evaluate_mesh_function(mesh, id2, boundary, qdim=2)

    if any(abs.(evaluate_function(vf, extract_nodes(boundary), qdim=2) .- vm) .> 1e-14)
        return false
    end

    for order = 1:8
        xle = length(quadrature_points_boundary(mesh.d, order))

        vq = evaluate_quadrature_function_boundary(mesh, id2, boundary, order=order, qdim=2)

        E = assemble_basismatrix_boundary(
            mesh,
            order = order,
            boundaryElements = belems, 
            qdim = 2
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
    belems = extract_elements(boundary)

    vf = evaluate_mesh_function(mesh, id3, qdim=3)
    vm = evaluate_mesh_function(mesh, id3, boundary, qdim=3)

    if any(abs.(evaluate_function(vf, extract_nodes(boundary), qdim=3) .- vm) .> 1e-14)
        return false
    end

    for order = 1:7
        points = quadrature_points_boundary(mesh, order)
        xle = length(quadrature_points_boundary(mesh.d, order))

        (
            points != quadrature_points_boundary(
                mesh,
                order,
                boundaryElements = belems
            )
        ) && return false

        vq = evaluate_quadrature_function_boundary(mesh, id3, boundary, order=order, qdim=3)

        E = assemble_basismatrix_boundary(
            mesh,
            order = order,
            boundaryElements = belems,
            qdim = 3
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
    belems = extract_elements(boundary)

    vm = evaluate_mesh_function(mesh, id3, boundary, qdim=3)

    if any(abs.(evaluate_function(vf, extract_nodes(boundary), qdim=3) .- vm) .> 1e-14)
        return false
    end

    for order = 1:7
        xle = length(quadrature_points_boundary(mesh.d, order))

        vq = evaluate_quadrature_function_boundary(mesh, id3, boundary, order=order, qdim=3)

        E = assemble_basismatrix_boundary(
            mesh,
            order = order,
            boundaryElements = belems, 
            qdim = 3
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
