using MinFEM

function test_output()
    mesh1::Mesh = import_mesh("test_line_v4.msh")
    mesh2::Mesh = import_mesh("test_square_v4.msh")
    mesh3::Mesh = import_mesh("test_cube_v4.msh")
    
    f11(x) = x[1]
    f21(x) = x[1] + x[2]
    f22(x) = [x[1]+x[2], x[1]-x[2]]
    f31(x) = x[1] + x[2] + x[3]
    f32(x) = [x[1]+x[2]+x[3], x[1]+x[2]-x[3]]
    f33(x) = [x[1]+x[2]+x[3], x[1]+x[2]-x[3], x[1]-x[2]-x[3]]

    fh11 = evaluate_mesh_function(mesh1, f11, qdim=1)
    fh21 = evaluate_mesh_function(mesh2, f21, qdim=1)
    fh22 = evaluate_mesh_function(mesh2, f22, qdim=2)
    fh31 = evaluate_mesh_function(mesh3, f31, qdim=1)
    fh32 = evaluate_mesh_function(mesh3, f32, qdim=2)
    fh33 = evaluate_mesh_function(mesh3, f33, qdim=3)

    pathBase = "temp"
    testDirPath = "" * pathBase
    while ispath(testDirPath)
        testDirPath = pathBase * "_$(rand(10000:100000))"
    end
    mkpath(testDirPath)

    write_to_txt(fh11, mesh1, testDirPath * "/test_f11"; qdim=1)
    write_to_txt(fh21, mesh2, testDirPath * "/test_f21"; qdim=1)
    write_to_txt(fh22, mesh2, testDirPath * "/test_f22.txt"; qdim=2)
    write_to_txt(fh31, mesh3, testDirPath * "/test_f31"; qdim=1)
    write_to_txt(fh32, mesh3, testDirPath * "/test_f32.txt"; qdim=2)
    write_to_txt(fh33, mesh3, testDirPath * "/test_f33.txt"; qdim=3)

    N11, v11 = read_from_txt(testDirPath * "/test_f11.txt")
    N21, v21 = read_from_txt(testDirPath * "/test_f21.txt")
    N22, v22 = read_from_txt(testDirPath * "/test_f22.txt")
    N31, v31 = read_from_txt(testDirPath * "/test_f31.txt")
    N32, v32 = read_from_txt(testDirPath * "/test_f32.txt")
    N33, v33 = read_from_txt(testDirPath * "/test_f33.txt")

    rm(testDirPath, recursive=true, force=true)

    for k = 1:mesh1.nnodes
        any(abs.(N11[k] .- mesh1.Nodes[k]) .> 1e-15) && return false
    end
    for k = 1:mesh2.nnodes
        any(abs.(N21[k] .- mesh2.Nodes[k]) .> 1e-15) && return false
        any(abs.(N22[k] .- mesh2.Nodes[k]) .> 1e-15) && return false
    end
    for k = 1:mesh3.nnodes
        any(abs.(N31[k] .- mesh3.Nodes[k]) .> 1e-15) && return false
        any(abs.(N32[k] .- mesh3.Nodes[k]) .> 1e-15) && return false
        any(abs.(N33[k] .- mesh3.Nodes[k]) .> 1e-15) && return false
    end

    any(abs.(v11 .- fh11) .> 1e-15) && return false
    any(abs.(v21 .- fh21) .> 1e-15) && return false
    any(abs.(v22 .- fh22) .> 1e-15) && return false
    any(abs.(v31 .- fh31) .> 1e-15) && return false
    any(abs.(v32 .- fh32) .> 1e-15) && return false
    any(abs.(v33 .- fh33) .> 1e-15) && return false

    return true
end

@test test_output()
