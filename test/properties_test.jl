using MinFEM

function test_properties_1d()
    mesh = unit_interval(3)
    bcs = elementbarycenter(mesh)
    any(abs.(bcs[1] .- 0.25) .> 1e-15) && return false
    any(abs.(bcs[2] .- 0.75) .> 1e-15) && return false
    for i = 1:mesh.nelems
        any(abs.(bcs[i] .- elementbarycenter(mesh, i)) .> 1e-15) && return false
    end
    any(x -> abs(x - 1.0) > 1e-15, elementratio(mesh)) && return false
    any(abs.(inscribedball(mesh) .- (1/4)) .> 1e-14) && return false
    any(abs.(circumscribedball(mesh) .- (1/4)) .> 1e-14) && return false

    mesh = unit_interval(11)
    vols = elementvolume(mesh)
    diams = elementdiameter(mesh)
    for i = 1:mesh.nelems
        abs(vols[i] - 1/10) > 1e-15 && return false
        abs(diams[i] - 1/10) > 1e-15 && return false
        abs(vols[i] - elementvolume(mesh, i)) > 1e-15 && return false
    end
    bvols = elementvolume_boundary(mesh)
    bdiams = elementdiameter_boundary(mesh)
    for i = 1:mesh.nboundelems
        abs(bvols[i] - 1) > 1e-15 && return false
        abs(bdiams[i] - 0) > 1e-15 && return false
        abs(bvols[i] - elementvolume_boundary(mesh, i)) > 1e-13 && return false
    end
    abs(gridsize(mesh) - 0.1) > 1e-14 && return false
    abs(shaperegularity(mesh) - 2.0) > 1e-14 && return false
    abs(quasiuniformity(mesh) - 1.0) > 1e-14 && return false

    mesh = import_mesh("test_line_v4.msh")
    abs(volume(mesh) - 1) > 1e-15 && return false
    any(abs.(barycenter(mesh) .- [0.5]) .> 1e-15) && return false
    bb = boundingbox(mesh) .- [[0], [1]]
    for k = 1:2
        for j = 1:mesh.d
            abs(bb[k][j]) > 1e-15 && return false
        end
    end
    abs(stripwidth(mesh) - 1) > 1e-15 && return false

    reg_interval = [[0.0],[1.0]]
    abs(inscribedball(reg_interval) - (1/2)) > 1e-14 && return false
    abs(circumscribedball(reg_interval) - (1/2)) > 1e-14 && return false
    abs(elementratio(reg_interval) - 1.0) > 1e-14 && return false

    return true
end

function test_properties_2d()
    mesh = unit_square(2)
    bcs = elementbarycenter(mesh)
    any(abs.(bcs[1] .- [1/3,1/3]) .> 1e-15) && return false
    any(abs.(bcs[2] .- [2/3,2/3]) .> 1e-15) && return false
    for i = 1:mesh.nelems
        any(abs.(bcs[i] .- elementbarycenter(mesh, i)) .> 1e-15) && return false
    end
    any(x -> abs(x - (sqrt(2)+1)) > 1e-14, elementratio(mesh)) && return false
    any(x -> abs(x - (45/180*pi)) > 1e-15, elementangle(mesh)) && return false

    mesh = unit_square(11)
    vols = elementvolume(mesh)
    diams = elementdiameter(mesh)
    for i = 1:mesh.nelems
        abs(vols[i] - 1/(2*10^2)) > 1e-15 && return false
        abs(diams[i] - sqrt(2)/10) > 1e-13 && return false
        abs(vols[i] - elementvolume(mesh, i)) > 1e-15 && return false
    end
    bvols = elementvolume_boundary(mesh)
    bdiams = elementdiameter_boundary(mesh)
    for i = 1:mesh.nboundelems
        abs(bvols[i] - 1/10) > 1e-15 && return false
        abs(bdiams[i] - 1/10) > 1e-15 && return false
        abs(bvols[i] - elementvolume_boundary(mesh, i)) > 1e-15 && return false
    end
 
    abs(gridsize(mesh) - 0.1*sqrt(2)) > 1e-14 && return false
    abs(shaperegularity(mesh) - (2+2*sqrt(2))) > 1e-14 && return false
    abs(quasiuniformity(mesh) - 1.0) > 1e-14 && return false

    mesh = import_mesh("test_square_v4.msh")
    abs(volume(mesh) - 1) > 1e-15 && return false
    any(abs.(barycenter(mesh) .- [0.5, 0.5]) .> 1e-15) && return false
    bb = boundingbox(mesh) .- [[0,0], [1,1]]
    for k = 1:2
        for j = 1:mesh.d
            abs(bb[k][j]) > 1e-15 && return false
        end
    end
    abs(stripwidth(mesh) - 1) > 1e-15 && return false

    reg_triangle = [[0.0,0.0],[1.0,0.0],[0.0,1.0]]
    abs(inscribedball(reg_triangle) - (1/(2+sqrt(2)))) > 1e-14 && return false
    abs(circumscribedball(reg_triangle) - (sqrt(2)/2)) > 1e-14 && return false
    abs(elementratio(reg_triangle) - (sqrt(2)+1)) > 1e-14 && return false

    return true
end

function test_properties_3d()
    mesh = import_mesh("test_cube_v4.msh")
    abs(volume(mesh) - 1) > 1e-14 && return false
    any(abs.(barycenter(mesh) .- [0.5, 0.5, 0.5]) .> 1e-14) && return false
    bb = boundingbox(mesh) .- [[0,0,0], [1,1,1]]
    for k = 1:2
        for j = 1:mesh.d
            abs(bb[k][j]) > 1e-15 && return false
        end
    end
    abs(stripwidth(mesh) - 1) > 1e-15 && return false
    any(
        x -> (
            abs(x - 45/180*pi) > 1e-15 &&
            abs(x - 60/180*pi) > 1e-15 &&
            abs(x - 30/180*pi) > 1e-15
        ),
        elementangle(mesh)
    ) && return false

    reg_tetrahedron = [[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]
    abs(inscribedball(reg_tetrahedron) - (1/(3+sqrt(3)))) > 1e-14 && return false
    abs(circumscribedball(reg_tetrahedron) - (sqrt(3)/2)) > 1e-14 && return false
    abs(elementratio(reg_tetrahedron) - (3/2*(sqrt(3)+1))) > 1e-14 && return false

    return true
end

function test_properties_fallbacks()
    try
        r = inscribedball([[0.0],[0.0],[0.0],[0.0],[0.0]])
    catch e
        if !isa(e, ArgumentError) || !occursin("Unsuitable set of coordinates", e.msg)
            return false
        end
    end
    try
        r = circumscribedball([[0.0],[0.0],[0.0],[0.0],[0.0]])
    catch e
        if !isa(e, ArgumentError) || !occursin("Unsuitable set of coordinates", e.msg)
            return false
        end
    end
    try
        r = elementratio([[0.0],[0.0],[0.0],[0.0],[0.0]])
    catch e
        if !isa(e, ArgumentError) || !occursin("Unsuitable set of coordinates", e.msg)
            return false
        end
    end
    
    try
        angles = elementangle(unit_interval(3))
    catch e
        if !isa(e, ErrorException) || !occursin("angle for 1D", e.msg)
            return false
        end
    end

    return true
end

@test test_properties_1d()
@test test_properties_2d()
@test test_properties_3d()
@test test_properties_fallbacks()