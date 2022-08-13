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

    export_mesh(meshi, testDirPath * "/testi.msh2")
    export_mesh(meshs, testDirPath * "/tests.msh2")
    remeshi = import_mesh(testDirPath * "/testi.msh2")
    remeshs = import_mesh(testDirPath * "/tests.msh2")

    rm(testDirPath, recursive=true, force=true)

    if !(remeshi == meshi && remeshs == meshs)
        return false
    end

    return true
end

@test test_mesh()
@test test_mesh("test_line")
@test test_mesh("test_square")
@test test_mesh("test_cube")
