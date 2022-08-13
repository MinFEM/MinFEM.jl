"""
    write_to_vtk(x::Vector{Vector{Float64}}, mesh::Mesh, dataName::Array{String}, fileName::String, qdim::Array{Int64})
    write_to_vtk(x::Vector{Vector{Float64}}, mesh::Mesh, dataNames::Array{String}, fileName::String; qdim::Int64=1)
    write_to_vtk(x::Vector{Float64}, mesh::Mesh, dataName::String, fileName::String; qdim=1)
    
Writes one or multiple vectors with variable qdim corrsponding to the mesh nodes to a VTK file. 
"""
function write_to_vtk(x::Vector{Vector{Float64}}, mesh::Mesh, dataName::Array{String}, fileName::String, qdim::Array{Int64})
    if !endswith(fileName,".vtu")
        fileName = fileName * ".vtu"
    end
    vtkfile = open_vtkfile(mesh, fileName)
    
    for k in eachindex(x)
        if qdim[k] == 1
            write_pointdata_vtkfile!(vtkfile, x[k], dataName[k])
        elseif qdim[k] == 2
            # Add third dimension to be able to use orientation array in paraview 
            write_pointdata_vtkfile!(vtkfile, [reshape(x[k], qdim[k], mesh.nnodes); zeros(1,mesh.nnodes)], dataName[k])
        else
            write_pointdata_vtkfile!(vtkfile, reshape(x[k], qdim[k], mesh.nnodes), dataName[k])
        end
    end
    vtk_save(vtkfile)
end

function write_to_vtk(x::Vector{Vector{Float64}}, mesh::Mesh, dataNames::Array{String}, fileName::String; qdim::Int64=1)
    write_to_vtk(x, mesh, dataNames, fileName, ones(Int64, length(x)) * qdim)
end

function write_to_vtk(x::Vector{Float64}, mesh::Mesh, dataName::String, fileName::String; qdim=1)
    write_to_vtk([x], mesh, [dataName], fileName, qdim=qdim)
end

"""
    write_to_vtk_boundary(x::Vector{Vector{Float64}}, mesh::Mesh, dataName::Array{String}, fileName::String, qdim::Array{Int64}; boundary=Set{Boundary}())
    write_to_vtk_boundary(x::Vector{Vector{Float64}}, mesh::Mesh, dataNames::Array{String}, fileName::String; boundary=Set{Boundary}(), qdim::Int64=1)
    write_to_vtk_boundary(x::Vector{Float64}, mesh::Mesh, dataName::String, fileName::String; boundary=Set{Boundary}(), qdim=1)
    
Writes one or multiple vectors with variable qdim corrsponding to the mesh boundary elements to a VTK file. 
"""
function write_to_vtk_boundary(x::Vector{Vector{Float64}}, mesh::Mesh, dataName::Array{String}, fileName::String, qdim::Array{Int64}; boundary=Set{Boundary}())
    if isempty(boundary)
        boundaryElements = Set{Int64}(1 : mesh.nboundelems)
    else
        boundaryElements = extract_elements(boundary)
    end
    
    if !endswith(fileName,".vtu")
        fileName = fileName * ".vtu"
    end
    vtkfile = open_vtkfile_boundary(mesh, fileName, boundaryElements)
    
    for k in eachindex(x)
        if qdim[k] == 1 
            val = x[k]
        elseif qdim[k] == 2
            # Add third dimension to be able to use orientation array in paraview
            val = [reshape(x[k], qdim[k], mesh.nboundelems); zeros(1,mesh.nboundelems)]
        else
            val = reshape(x[k], qdim[k], mesh.nboundelems)
        end
        val = val[:,1:end .∈ [boundaryElements]]

        write_celldata_vtkfile!(vtkfile, val, dataName[k])
    end
    vtk_save(vtkfile)
end

function write_to_vtk_boundary(x::Vector{Vector{Float64}}, mesh::Mesh, dataNames::Array{String}, fileName::String; boundary=Set{Boundary}(), qdim::Int64=1)
    write_to_vtk_boundary(x, mesh, dataNames, fileName, ones(Int64, length(x)) * qdim, boundary=boundary)
end

function write_to_vtk_boundary(x::Vector{Float64}, mesh::Mesh, dataName::String, fileName::String; boundary=Set{Boundary}(), qdim=1)
    write_to_vtk_boundary([x], mesh, [dataName], fileName, qdim=qdim, boundary=boundary)
end

"""
    open_vtkfile(mesh::Mesh, file_name::String) -> WriteVTK.DatasetFile

Open a new VTK output file and write the mesh data into it.
"""
function open_vtkfile(mesh::Mesh, file_name::String)
    points = zeros(Float64, length(mesh.Nodes[1]), length(mesh.Nodes))
    for (i,p) in enumerate(mesh.Nodes)
        points[:,i] = copy(p)
    end

    cells = Array{MeshCell,1}(undef, mesh.nelems)
    if mesh.d == 1
        for (i,t) in enumerate(mesh.Elements)
            cells[i] = MeshCell(VTKCellTypes.VTK_LINE, t)
        end   
    elseif mesh.d == 2
        for (i,t) in enumerate(mesh.Elements)
            cells[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, t)
        end
    elseif mesh.d == 3
        for (i,t) in enumerate(mesh.Elements)
            cells[i] = MeshCell(VTKCellTypes.VTK_TETRA, t)
        end
    end

    return vtk_grid(file_name, points, cells)
end

"""
    open_vtkfile_boundary(mesh::Mesh, file_name::String, boundaryElements=Set{Int64}()) -> WriteVTK.DatasetFile
    open_vtkfile_boundary(mesh::Mesh, file_name::String; boundary=Set{Boundary}()) -> WriteVTK.DatasetFile

Open a new VTK output file and write the mesh data into it.
"""
function open_vtkfile_boundary(mesh::Mesh, file_name::String, boundaryElements::Set{Int64})
    if isempty(boundaryElements)
        boundaryElements = Set{Int64}(1 : mesh.nboundelems)
    end
    support = mesh.BoundaryElements[1:end .∈ [boundaryElements]]

    points = zeros(Float64, length(mesh.Nodes[1]), length(mesh.Nodes))
    for (i,p) in enumerate(mesh.Nodes)
        points[:,i] = copy(p)
    end

    if mesh.d == 1
        type = VTKCellTypes.VTK_VERTEX
    elseif mesh.d == 2
        type = VTKCellTypes.VTK_LINE
    elseif mesh.d == 3
        type = VTKCellTypes.VTK_TRIANGLE
    else
        throw(DomainError(mesh.d, "Invalid dimension."))
    end

    cells = Array{MeshCell,1}(undef, length(support))
    for (i,t) in enumerate(support)
        cells[i] = MeshCell(type, t)
    end

    return vtk_grid(file_name, points, cells)
end

function open_vtkfile_boundary(mesh::Mesh, file_name::String; boundary=Set{Boundary}())
    open_vtkfile_boundary(mesh, file_name, extract_elements(boundary))
end

"""
    write_pointdata_vtkfile!(vtkfile::WriteVTK.DatasetFile, data::Any, data_name::String)

Add a new point data field with a name to an existing VTK file.
"""
function write_pointdata_vtkfile!(vtkfile::WriteVTK.DatasetFile, data::Any, data_name::String)
    WriteVTK.vtk_point_data(vtkfile, data, data_name)
end

"""
    write_celldata_vtkfile!(vtkfile::WriteVTK.DatasetFile, data::Any, data_name::String)  

Add a new cell data field with a name to an existing VTK file.
"""
function write_celldata_vtkfile!(vtkfile::WriteVTK.DatasetFile, data::Any, data_name::String)
    WriteVTK.vtk_cell_data(vtkfile, data, data_name)
end

"""
    save_vtkfile(vtkfile::WriteVTK.DatasetFile)

Finalize a VTK file by writing all data to disk.
"""
function save_vtkfile(vtkfile::WriteVTK.DatasetFile)
    WriteVTK.vtk_save(vtkfile)
end


"""
    write_to_txt(x::Vector{Float64}, mesh::Mesh, fileName::String; qdim::Int64=1)

Writes a coefficient vector x based on the nodes of mesh to the given file.
"""
function write_to_txt(x::Vector{Float64}, mesh::Mesh, fileName::String; qdim::Int64=1)
    if !endswith(fileName,".txt")
        fileName = fileName * ".txt"
    end

    open(fileName, "w") do file
        write(file,"\$PARAMS\n")
        write(file,"$(mesh.d) $qdim $(mesh.nnodes)\n")
        write(file, "\$VALS\n")
        for i = 1:mesh.nnodes
            write(file, string(i))
            for j = 1:mesh.d
                write(file, " ", string(mesh.Nodes[i][j]))
            end

            for r = 1:qdim
                write(file, " ", string(x[qdim*(i-1)+r]))
            end
            write(file,"\n")
        end
    end
end

"""
    read_from_txt(fileName::String)

Reads node coordinates and finite element coefficient vector from file generated by `write_to_txt`.
"""
function read_from_txt(fileName::String)
    f = open(fileName)
    while(!eof(f) && (l=readline(f)) != "\$PARAMS") end
    
    l = readline(f)
    a = split(l, " ")
    dim = parse(Int64, a[1])
    qdim = parse(Int64, a[2])
    nvals = parse(Int64, a[3])

    Nodes = Array{Array{Float64},1}(undef, nvals)
    v = Array{Float64, 1}(undef, nvals*qdim)

    while(!eof(f) && (l=readline(f)) != "\$VALS") end
    for i = 1:nvals
        l = readline(f)
        a = split(l, " ")

        coords = Array{Float64,1}(undef, dim)
        for j = 1:dim
            coords[j] = parse(Float64, a[j+1])
        end
        Nodes[i] = coords

        for r = 1:qdim
            v[qdim*(i-1)+r] = parse(Float64, a[1+dim+r])
        end
    end
    
    close(f)

    return Nodes, v
end
