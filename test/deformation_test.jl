using MinFEM

function test_deformation()
    f1(x) = max(x[1]-0.5, 0)
    f2(x) = [x[1]-0.5, x[2]-0.5]
    f3(x) = [1, 0, -1] 
    
    mesh = unit_interval(11)
    fh = evaluate_mesh_function(mesh, f1)

    mesh2 = deform_mesh(mesh, fh, t=2.0)
    deform_mesh!(mesh, fh, t=2.0)
    mesh != mesh2 && return false

    mesh3 = unit_interval(11)
    update_mesh!(mesh3, mesh2.Nodes)
    mesh != mesh3 && return false

    diff = mesh.Nodes .- [[0.0], [0.1], [0.2], [0.3], [0.4], [0.5], [0.8], [1.1], [1.4], [1.7], [2.0]]
    for k = 1:mesh.nnodes
        abs(diff[k][1]) > 1e-15 && return false
    end
    abs(volume(mesh) - 2) > 1e-15 && return false
    any(abs.(barycenter(mesh) .- [1]) .> 1e-15) && return false

    
    mesh = unit_square(11)
    fh = evaluate_mesh_function(mesh, f2, qdim=2)

    mesh2 = deform_mesh(mesh, fh, t=2.0)
    deform_mesh!(mesh, fh, t=2.0)
    mesh != mesh2 && return false

    mesh3 = unit_square(11)
    update_mesh!(mesh3, mesh2.Nodes)
    mesh != mesh3 && return false

    abs(volume(mesh) - 9) > 1e-13 && return false
    any(abs.(barycenter(mesh) .- [0.5, 0.5]) .> 1e-15) && return false
    bb = boundingbox(mesh) .- [[-1,-1], [2,2]]
    for k = 1:2
        any(abs.(bb[k]) .> 1e-15) && return false
    end


    mesh = import_mesh("test_cube_v4.msh")
    fh = evaluate_mesh_function(mesh, f3, qdim=3)

    mesh2 = deform_mesh(mesh, fh)
    deform_mesh!(mesh, fh)
    mesh != mesh2 && return false

    mesh3 = import_mesh("test_cube_v4.msh")
    update_mesh!(mesh3, mesh2.Nodes)
    mesh != mesh3 && return false
    
    abs(volume(mesh) - 1) > 1e-13 && return false
    any(abs.(barycenter(mesh) .- [1.5, 0.5, -0.5]) .> 1e-14) && return false
    bb = boundingbox(mesh) .- [[1,0,-1], [2,1,0]]
    for k = 1:2
        any(abs.(bb[k]) .> 1e-15) && return false
    end

    return true
end

function test_deformation_fallbacks()
    mesh = unit_interval(11)
    fh = zeros(Float64,mesh.nelems)
    try
        deform_mesh!(mesh, fh)
    catch e
        if !isa(e, ArgumentError) || !occursin("matching length", e.msg)
            return false
        end
    end

    mesh = unit_square(11)
    fh = zeros(Float64,mesh.nnodes)
    try
        deform_mesh!(mesh, fh)
    catch e
        if !isa(e, ArgumentError) || !occursin("matching length", e.msg)
            return false
        end
    end

    mesh = import_mesh("test_cube_v4.msh")
    fh = zeros(Float64, mesh.nnodes * 2)
    try
        deform_mesh!(mesh, fh)
    catch e
        if !isa(e, ArgumentError) || !occursin("matching length", e.msg)
            return false
        end
    end

    return true
end

@test test_deformation()
@test test_deformation_fallbacks()
