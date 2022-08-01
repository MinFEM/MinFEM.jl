using MinFEMDev

function test_properties()
    mesh = unit_interval(11)
    vols = elementvolume(mesh)
    diams = elementdiameter(mesh)
    for i = 1:mesh.nelems
        abs(vols[i] - 1/10) > 1e-15 && return false
        abs(diams[i] - 1/10) > 1e-15 && return false
    end
    bvols = elementvolume_boundary(mesh)
    bdiams = elementdiameter_boundary(mesh)
    for i = 1:mesh.nboundelems
        abs(bvols[i] - 1) > 1e-15 && return false
        abs(bdiams[i] - 0) > 1e-15 && return false
    end

    mesh = unit_square(11)
    vols = elementvolume(mesh)
    diams = elementdiameter(mesh)
    for i = 1:mesh.nelems
        abs(vols[i] - 1/(2*10^2)) > 1e-15 && return false
        abs(diams[i] - sqrt(2)/10) > 1e-13 && return false
    end
    bvols = elementvolume_boundary(mesh)
    bdiams = elementdiameter_boundary(mesh)
    for i = 1:mesh.nboundelems
        abs(bvols[i] - 1/10) > 1e-15 && return false
        abs(bdiams[i] - 1/10) > 1e-15 && return false
    end

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

    return true
end

@test test_properties()
