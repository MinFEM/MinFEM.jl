using MinFEM

function test_mesh(filePath::String)
    mesh1 = import_mesh("$(filePath)_v1.msh")
    mesh2 = import_mesh("$(filePath)_v2.msh")
    mesh4 = import_mesh("$(filePath)_v4.msh")

    if !(mesh1 == mesh2 && mesh2 == mesh4)
        return false
    end

    pathBase = "temp"
    testDirPath = "" * pathBase
    while ispath(testDirPath)
        testDirPath = pathBase * "_$(rand(10000:100000))"
    end
    mkpath(testDirPath)

    export_mesh(mesh1, testDirPath * "/test_v1.msh")
    export_mesh(mesh2, testDirPath * "/test_v2.msh")
    export_mesh(mesh4, testDirPath * "/test_v4.msh")
    remesh1 = import_mesh(testDirPath * "/test_v1.msh")
    remesh2 = import_mesh(testDirPath * "/test_v2.msh")
    remesh4 = import_mesh(testDirPath * "/test_v4.msh")

    rm(testDirPath, recursive=true, force=true)

    if !(remesh1 == remesh2 && remesh2 == remesh4)
        return false
    end

    if !(remesh2 == mesh2)
        return false
    end

    return true

end

function test_mesh()
    meshi = unit_interval(10)
    meshs = unit_square(4)

    pathBase = "temp"
    testDirPath = "" * pathBase
    while ispath(testDirPath)
        testDirPath = pathBase * "_$(rand(10000:100000))"
    end
    mkpath(testDirPath)

    export_mesh(meshi, testDirPath * "/testi")
    export_mesh(meshs, testDirPath * "/tests")
    remeshi = import_mesh(testDirPath * "/testi.msh2")
    remeshs = import_mesh(testDirPath * "/tests.msh2")

    rm(testDirPath, recursive=true, force=true)

    if !(remeshi == meshi && remeshs == meshs)
        return false
    end

    return true
end

function test_mesh_fallbacks()
    if unit_interval(1).nnodes != 2
        return false
    end

    if unit_square(1).nnodes != 4
        return false
    end

    try
        import_mesh("test_error_empty.msh")
    catch e
        if !isa(e, ErrorException) || !occursin("corrupted", e.msg)
            return false
        end
    end

    try
        import_mesh("test_error_v3.msh")
    catch e
        if !isa(e, ErrorException) || !occursin("Unsupported", e.msg)
            return false
        end
    end

    try
        import_mesh("test_error_v2.msh")
    catch e
        if !isa(e, ErrorException) || !occursin("meshed domain", e.msg)
            return false
        end
    end

    try
        import_mesh("test_error_v4.msh")
    catch e
        if !isa(e, ErrorException) || !occursin("meshed domain", e.msg)
            return false
        end
    end

    return true
end

function test_mesh_jacobians()
    mesh = import_mesh("test_line_v4.msh")
    el = 1
    j1,_ = jacobian(mesh, el)
    j2,_ = jacobian(mesh, mesh.Elements[el])
    j3,_ = jacobian(mesh.Nodes[mesh.Elements[el]])
    (j1 != j2) && return false
    (j2 != j3) && return false

    mesh = import_mesh("test_square_v4.msh")
    el = mesh.nelems
    j1,_ = jacobian(mesh, el)
    j2,_ = jacobian(mesh, mesh.Elements[el])
    j3,_ = jacobian(mesh.Nodes[mesh.Elements[el]])
    (j1 != j2) && return false
    (j2 != j3) && return false

    mesh = import_mesh("test_cube_v4.msh")
    el = 1
    j1,_ = jacobian(mesh, el)
    j2,_ = jacobian(mesh, mesh.Elements[el])
    j3,_ = jacobian(mesh.Nodes[mesh.Elements[el]])
    (j1 != j2) && return false
    (j2 != j3) && return false
    
    return true
end

function test_mesh_regions()
    mesh = unit_square(3)

    ball = select_boundaries(mesh)
    (ball != select_boundaries(mesh, 1001, 1002, 1003, 1004)) && return false

    b1001 = select_boundaries(mesh, 1001)
    try
        b1005 = select_boundaries(mesh, 1005)
    catch e
        if !isa(e, ArgumentError) || !occursin("boundary id", e.msg)
            return false
        end
    end

    !isa(b1001, Set{Boundary}) && return false
    (extract_elements(b1001) != Set{Int64}(1:2)) && return false
    (extract_nodes(b1001) != Set{Int64}(1:3)) && return false

    d10001 = select_domains(mesh, 10001)
    try
        d10002 = select_domains(mesh, 10002)
    catch e
        if !isa(e, ArgumentError) || !occursin("domain id", e.msg)
            return false
        end
    end

    !isa(d10001, Set{Domain}) && return false
    (d10001 != select_domains(mesh)) && return false
    (extract_elements(d10001) != Set{Int64}(1:8)) && return false
    (extract_nodes(d10001) != Set{Int64}(1:9)) && return false    
    
    return true
end

function test_mesh_quality()
    mesh = unit_interval(3)
    ratios = elementratio(mesh)
    any(x -> abs(x - 1.0) > 1e-5, ratios) && return false
    try
        angles = elementangle(mesh)
    catch e
        if !isa(e, ErrorException) || !occursin("angle for 1D", e.msg)
            return false
        end
    end
    
    mesh = unit_square(3)
    ratios = elementratio(mesh)
    any(x -> abs(x - 0.41421) > 1e-5, ratios) && return false
    angles = elementangle(mesh)
    any(x -> abs(x - 0.78539) > 1e-5, angles) && return false

    mesh = import_mesh("test_cube_v4.msh")
    try
        ratios = elementratio(mesh)
    catch e
        if !isa(e, ErrorException) || !occursin("ratio for 3D", e.msg)
            return false
        end
    end
    angles = elementangle(mesh)
    any(
        x -> (
            abs(x - 0.78539) > 1e-5 &&
            abs(x - 1.04719) > 1e-5 &&
            abs(x - 0.52359) > 1e-5
        ),
        angles
    ) && return false

    return true
end

@test test_mesh()
@test test_mesh("test_line")
@test test_mesh("test_square")
@test test_mesh("test_cube")
@test test_mesh_fallbacks()
@test test_mesh_jacobians()
@test test_mesh_regions()
@test test_mesh_quality()
