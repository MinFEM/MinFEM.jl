using MinFEMDev

function test_normals()
    f(x) = 0.5 * x[1]
    g(x) = 0.5 .* [x[1], -x[3], x[2]]

    mesh::Mesh = import_mesh("test_line_v4.msh")
    belems = extract_elements(select_boundaries(mesh))
    fh = evaluate_mesh_function(mesh, f)

    n = outernormalvector(mesh, boundaryElements=belems)
    any(abs.(n .- [-1.0, 1.0]) .> 1e-15) && return false

    Dn = assemble_normalderivativematrix(mesh, boundaryElements=belems, qdim=1)
    Dnf = Dn * fh
    any(abs.(Dnf .- [-0.5, 0.5]) .> 1e-15) && return false


    mesh = import_mesh("test_square_v4.msh")
    belems = extract_elements(select_boundaries(mesh))
    fh = evaluate_mesh_function(mesh, f)

    soln = zeros(mesh.nboundelems*2)
    solDnf = zeros(mesh.nboundelems)
    for bel in extract_elements(select_boundaries(mesh, 1001))
        soln[2*(bel-1)+1 : 2*bel] = [0,-1]
    end
    for bel in extract_elements(select_boundaries(mesh, 1002))
        soln[2*(bel-1)+1 : 2*bel] = [0,1]
    end
    for bel in extract_elements(select_boundaries(mesh, 1003))
        soln[2*(bel-1)+1 : 2*bel] = [-1,0]
        solDnf[bel] = -0.5
    end
    for bel in extract_elements(select_boundaries(mesh, 1004))
        soln[2*(bel-1)+1 : 2*bel] = [1,0]
        solDnf[bel] = 0.5
    end

    n = outernormalvector(mesh, boundaryElements=belems)
    any(abs.(n .- soln) .> 1e-15) && return false

    Dn = assemble_normalderivativematrix(mesh, boundaryElements=belems, qdim=1)
    Dnf = Dn * fh
    any(abs.(Dnf .- solDnf) .> 1e-15) && return false

    mesh = import_mesh("test_cube_v4.msh")
    belems = extract_elements(select_boundaries(mesh))
    fh = evaluate_mesh_function(mesh, f)
    gh = evaluate_mesh_function(mesh, g, qdim=3)

    soln = zeros(mesh.nboundelems*3)
    solDnf = zeros(mesh.nboundelems)
    solDng = zeros(mesh.nboundelems*3)
    for bel in extract_elements(select_boundaries(mesh, 1001))
        soln[3*(bel-1)+1 : 3*bel] = [0,0,1]
        solDng[3*(bel-1)+1 : 3*bel] = [0,-0.5,0]
    end
    for bel in extract_elements(select_boundaries(mesh, 1002))
        soln[3*(bel-1)+1 : 3*bel] = [1,0,0]
        solDnf[bel] = 0.5
        solDng[3*(bel-1)+1 : 3*bel] = [0.5,0,0]
    end
    for bel in extract_elements(select_boundaries(mesh, 1003))
        soln[3*(bel-1)+1 : 3*bel] = [0,0,-1]
        solDng[3*(bel-1)+1 : 3*bel] = [0,0.5,0]
    end
    for bel in extract_elements(select_boundaries(mesh, 1004))
        soln[3*(bel-1)+1 : 3*bel] = [-1,0,0]
        solDnf[bel] = -0.5
        solDng[3*(bel-1)+1 : 3*bel] = [-0.5,0,0]
    end
    for bel in extract_elements(select_boundaries(mesh, 1005))
        soln[3*(bel-1)+1 : 3*bel] = [0,-1,0]
        solDng[3*(bel-1)+1 : 3*bel] = [0,0,-0.5]
    end
    for bel in extract_elements(select_boundaries(mesh, 1006))
        soln[3*(bel-1)+1 : 3*bel] = [0,1,0]
        solDng[3*(bel-1)+1 : 3*bel] = [0,0,0.5]
    end

    n = outernormalvector(mesh, boundaryElements=belems)
    any(abs.(n .- soln) .> 1e-15) && return false

    Dn = assemble_normalderivativematrix(mesh, boundaryElements=belems)
    Dnf = Dn * fh
    any(abs.(Dnf .- solDnf) .> 1e-15) && return false

    Dn = assemble_normalderivativematrix(mesh, boundaryElements=belems, qdim=3)
    Dng = Dn * gh
    any(abs.(Dng .- solDng) .> 1e-15) && return false

    return true
end

@test test_normals()
