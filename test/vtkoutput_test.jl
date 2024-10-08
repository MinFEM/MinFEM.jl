using MinFEM

function test_vtk()
    mesh1::Mesh = import_mesh("test_line_v4.msh")
    mesh2::Mesh = import_mesh("test_square_v4.msh")
    mesh3::Mesh = import_mesh("test_cube_v4.msh")

    b1 = extract_elements(select_boundaries(mesh1))
    b2 = extract_elements(select_boundaries(mesh2))
    b3 = extract_elements(select_boundaries(mesh3))
    
    EV11 = assemble_basismatrix(mesh1, order=1, qdim=1)
    EV21 = assemble_basismatrix(mesh2, order=1, qdim=1)
    EV22 = assemble_basismatrix(mesh2, order=1, qdim=2)
    EV31 = assemble_basismatrix(mesh3, order=1, qdim=1)
    EV32 = assemble_basismatrix(mesh3, order=1, qdim=2)
    EV33 = assemble_basismatrix(mesh3, order=1, qdim=3)
    
    E11 = assemble_basismatrix_boundary(mesh1, boundaryElements=b1, order=1, qdim=1)
    E21 = assemble_basismatrix_boundary(mesh2, boundaryElements=b2, order=1, qdim=1)
    E22 = assemble_basismatrix_boundary(mesh2, boundaryElements=b2, order=1, qdim=2)
    E31 = assemble_basismatrix_boundary(mesh3, boundaryElements=b3, order=1, qdim=1)
    E32 = assemble_basismatrix_boundary(mesh3, boundaryElements=b3, order=1, qdim=2)
    E33 = assemble_basismatrix_boundary(mesh3, boundaryElements=b3, order=1, qdim=3)

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

    vfh11 = EV11 * fh11
    vfh21 = EV21 * fh21
    vfh22 = EV22 * fh22
    vfh31 = EV31 * fh31
    vfh32 = EV32 * fh32
    vfh33 = EV33 * fh33

    bfh11 = E11 * fh11
    bfh21 = E21 * fh21
    bfh22 = E22 * fh22
    bfh31 = E31 * fh31
    bfh32 = E32 * fh32
    bfh33 = E33 * fh33

    pathBase = "temp"
    testDirPath = "" * pathBase
    while ispath(testDirPath)
        testDirPath = pathBase * "_$(rand(10000:100000))"
    end
    mkpath(testDirPath)

    success::Bool = true

    try
        write_to_vtk(fh11, mesh1, "f", testDirPath * "/test_f11a")
        write_to_vtk(fh11, mesh1, "f", testDirPath * "/test_f11b", qdim=1)

        write_to_vtk(vfh11, mesh1, "f", testDirPath * "/test_vf11a")
        write_to_vtk(vfh11, mesh1, "f", testDirPath * "/test_vf11b", qdim=1)

        v1 = open_vtkfile_boundary(mesh1, testDirPath * "/test_bv1")
        save_vtkfile(v1)

        write_to_vtk_boundary(bfh11, mesh1, "f", testDirPath * "/test_bf11a")
        write_to_vtk_boundary(bfh11, mesh1, "f", testDirPath * "/test_bf11b", qdim=1)
        write_to_vtk_boundary(bfh11, mesh1, "f", testDirPath * "/test_bf11c", qdim=1, boundary=select_boundaries(mesh1))

        write_to_vtk(fh21, mesh2, "f", testDirPath * "/test_f21")
        write_to_vtk(fh22, mesh2, "f", testDirPath * "/test_f22a", qdim=2)

        write_to_vtk([fh22, fh22], mesh2, ["f22a","f22b"], testDirPath * "/test_f22b", qdim=2)
        write_to_vtk([fh22, fh21], mesh2, ["f22","f21"], testDirPath * "/test_f2", [2,1])

        write_to_vtk(vfh21, mesh2, "f", testDirPath * "/test_vf21")
        write_to_vtk(vfh22, mesh2, "f", testDirPath * "/test_vf22a", qdim=2)
        write_to_vtk([vfh22, vfh22], mesh2, ["f22a","f22b"], testDirPath * "/test_vf22b", qdim=2)
        write_to_vtk([vfh22, vfh21], mesh2, ["f22","f21"], testDirPath * "/test_vf2", [2,1])

        write_to_vtk_boundary(bfh21, mesh2, "f", testDirPath * "/test_bf21")
        write_to_vtk_boundary(bfh22, mesh2, "f", testDirPath * "/test_bf22a", qdim=2)
        write_to_vtk_boundary([bfh22, bfh22], mesh2, ["f22a","f22b"], testDirPath * "/test_bf22b", qdim=2)
        write_to_vtk_boundary([bfh22, bfh21], mesh2, ["f22","f21"], testDirPath * "/test_bf2", [2,1])

        write_to_vtk(fh31, mesh3, "f", testDirPath * "/test_f31")
        write_to_vtk(fh32, mesh3, "f", testDirPath * "/test_f32", qdim=2)
        write_to_vtk(fh33, mesh3, "f", testDirPath * "/test_f33a", qdim=3)
        write_to_vtk([fh33, fh33, fh33], mesh3, ["f33a","f33b","f33c"], testDirPath * "/test_f33b", qdim=3)
        write_to_vtk([fh33, fh32, fh31], mesh3, ["f33","f32","f31"], testDirPath * "/test_f3", [3,2,1])

        write_to_vtk(vfh31, mesh3, "f", testDirPath * "/test_vf31")
        write_to_vtk(vfh32, mesh3, "f", testDirPath * "/test_vf32", qdim=2)
        write_to_vtk(vfh33, mesh3, "f", testDirPath * "/test_vf33a", qdim=3)
        write_to_vtk([vfh33, vfh33, vfh33], mesh3, ["f33a","f33b","f33c"], testDirPath * "/test_vf33b", qdim=3)
        write_to_vtk([vfh33, vfh32, vfh31], mesh3, ["f33","f32","f31"], testDirPath * "/test_vf3", [3,2,1])

        write_to_vtk_boundary(bfh31, mesh3, "f", testDirPath * "/test_bf31")
        write_to_vtk_boundary(bfh32, mesh3, "f", testDirPath * "/test_bf32", qdim=2)
        write_to_vtk_boundary(bfh33, mesh3, "f", testDirPath * "/test_bf33a", qdim=3)
        write_to_vtk_boundary([bfh33, bfh33, bfh33], mesh3, ["f33a","f33b","f33c"], testDirPath * "/test_bf33b", qdim=3)
        write_to_vtk_boundary([bfh33, bfh32, bfh31], mesh3, ["f33","f32","f31"], testDirPath * "/test_bf3", [3,2,1])
    catch
        success = false
    end

    try
        fh11mod = fh11[1:(end-2)]
        write_to_vtk(fh11mod, mesh1, "f", testDirPath * "/test_f11mod")
        success = false # fail if error was not thrown
    catch e
        if !isa(e, DimensionMismatch) || !occursin("vector does not match", e.msg)
            success = false
        end
    end

    rm(testDirPath, recursive=true, force=true)

    return success
end

@test test_vtk()
