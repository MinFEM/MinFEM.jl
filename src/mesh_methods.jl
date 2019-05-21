"""
    unit_square(n::Int64)

Set up a n-by-n quasi-uniform mesh for the unit square.
"""
function unit_square(n::Int64)
  if(n<2)
    n = 2
  end
  nnodes = n*n
  nedges = 4*(n-1)
  nelems = 2*(n-1)*(n-1)
  h = 1/(n-1)

  Nodes = Array{Array{Float64, 1}, 1}(undef, nnodes)
  Edges = Array{Array{Int64,1},1}(undef, nedges)
  Triangles = Array{Array{Int64,1},1}(undef, nelems)
  Volumes = Dict{Int64,Set{Int64}}()
  Boundaries = Dict{Int64, Boundary}()
  Boundaries[1001] = Boundary(Set{Int64}(), Set{Int64}())
  Boundaries[1002] = Boundary(Set{Int64}(), Set{Int64}())
  Boundaries[1003] = Boundary(Set{Int64}(), Set{Int64}())
  Boundaries[1004] = Boundary(Set{Int64}(), Set{Int64}())

  k = 1
  for i=1:n
    for j=1:n
      Nodes[k] = [(j-1)*h; (i-1)*h]
      k = k+1
    end
  end

  k = 1
  for i=1:n-1
    for j=1:n-1
      Triangles[k] = [(i-1)*n + j; (i-1)*n + j + 1; i*n + j]
      k = k+1

      Triangles[k] = [(i-1)*n + j + 1; i*n + j + 1; i*n + j]
      k = k+1
    end
  end

  for j=1:n
    push!(Boundaries[1001].Nodes, j)
    push!(Boundaries[1002].Nodes, j + (n-1)*n)
    push!(Boundaries[1003].Nodes, (j-1)*n+1)
    push!(Boundaries[1004].Nodes, (j-1)*n+n)
  end

  k = 1
  for j=1:n-1
    Edges[k] =  [j, j+1];
    push!(Boundaries[1001].Edges, k)
    k = k+1

    Edges[k] =  [j+1 + (n-1)*n, j + (n-1)*n];
    push!(Boundaries[1002].Edges, k)
    k = k+1

    Edges[k] =  [j*n+1, (j-1)*n+1];
    push!(Boundaries[1003].Edges, k)
    k = k+1

    Edges[k] =  [(j-1)*n+n, j*n+n];
    push!(Boundaries[1004].Edges, k)
    k = k+1
  end

  Volumes[10001] = Set{Int64}(1:nnodes)

  return Mesh(nnodes, nelems, nedges, Nodes, Triangles, Edges, Boundaries, Volumes)
end

function import_mesh(file_name::String)
  f = open(file_name)

  while(!eof(f) && (l=readline(f)) != "\$MeshFormat")
  end
  l=readline(f)
  a = split(l, " ")
  version = parse(Float64, a[1])
  if(version >= 2.0 && version <3.0)
    mesh = import_mesh2(f)
  elseif(version >= 4.1 && version <5.0)
    mesh = import_mesh4(f)
  else
    println("Unsupported mesh format: " * string(parse(Float64, a[1])))
    println("msh2 is recommended.")
    return
  end
  close(f)
  return mesh
end

function import_mesh2(f::IOStream)

  while(!eof(f) && (l=readline(f)) != "\$Nodes")
  end
  l=readline(f)
  nnodes = parse(Int64, l)
  Nodes = zeros(nnodes)
  Nodes = Array{Array{Float64, 1}, 1}(undef, nnodes)

  for i=1:nnodes
    l = readline(f)
    a = split(l, " ")
    Nodes[i] = [parse(Float64, a[2]), parse(Float64, a[3])]
  end
  while(!eof(f) && (l=readline(f)) != "\$Elements")
  end
  l=readline(f)
  nelems = parse(Int64, l)

  Triangles = Set(Array{Int64}[])
  Edges = Set(Array{Int64}[])
  for i=1:nelems
    l = readline(f)
    a = split(l, " ")
    if(parse(Int64, a[2]) == 1)
      push!(Edges, [parse(Int64, a[4]), parse(Int64, a[6]),
                    parse(Int64, a[7])])
    else
      push!(Triangles, [parse(Int64, a[4]), parse(Int64, a[6]),
                        parse(Int64, a[7]), parse(Int64, a[8])])
    end
    
  end

  _Triangles = Array{Array{Int64,1},1}(undef,length(Triangles))
  _Edges = Array{Array{Int64,1},1}(undef,length(Edges))

  Volumes = Dict{Int64,Set{Int64}}()
  for (i,t) in enumerate(Triangles)
    _Triangles[i] = copy(t[2:end])
    if( !in(t[1],keys(Volumes)) )
      Volumes[t[1]] = Set{Int64}()
    end
    push!(Volumes[t[1]], i)
  end

  Boundaries = Dict{Int64, Boundary}()
  for (i,el) in enumerate(Edges)
    _Edges[i] = copy(el[2:end])
    if( !in(el[1],keys(Boundaries)) )
      Boundaries[el[1]] = Boundary(Set{Int64}(), Set{Int64}())
    end
    push!(Boundaries[el[1]].Edges, i)
    push!(Boundaries[el[1]].Nodes, el[2])
    push!(Boundaries[el[1]].Nodes, el[3])
  end

  return Mesh(nnodes, length(Triangles), length(Edges), Nodes, _Triangles, _Edges, Boundaries, Volumes)
end

function import_mesh4(f::IOStream)

  # 0,1,2,3 dimensional entity tags
  tags = [Dict{Int64,Int64}(), Dict{Int64,Int64}(),
          Dict{Int64,Int64}(), Dict{Int64,Int64}()]
  while(!eof(f) && (l=readline(f)) != "\$Entities")
  end
  l=readline(f)
  a=split(l, " ")
  numTags = [parse(Int64, a[1]), parse(Int64, a[2]),
             parse(Int64, a[3]), parse(Int64, a[4])]

  # Point tags
  for j=1:numTags[1]
    l=readline(f)
    a=split(l, " ")
    if(parse(Int64, a[4])!=0)
      tags[1][parse(Int64, a[1])] = parse(Int64, a[5])
    end
  end

  # Curve, Surface, Volume tags
  for i=2:4
    for j=1:numTags[i]
      l=readline(f)
      a=split(l, " ")
      if(parse(Int64, a[8])!=0)
        tags[i][parse(Int64, a[1])] = parse(Int64, a[9])
      end
    end
  end

  while(!eof(f) && (l=readline(f)) != "\$Nodes")
  end
  l=readline(f)
  a=split(l, " ")
  blocks = parse(Int64, a[1])
  nnodes = parse(Int64, a[2])

  Nodes = Array{Array{Float64, 1}, 1}(undef, nnodes)
  NodeNumbering = Dict{Int64, Int64}()
  n=1
  m=1
  for i=1:blocks
    l = readline(f)
    a = split(l, " ")
    nodesInBlock = parse(Int64, a[4])

    # The node numbers
    for j=1:nodesInBlock
      l = readline(f)
      a = split(l, " ")
      NodeNumbering[parse(Int64, a[1])] = n
      n+=1
    end

    # The actual coordinates
    for j=1:nodesInBlock
      l = readline(f)
      a = split(l, " ")
      Nodes[m] = [parse(Float64, a[1]), parse(Float64, a[2])]
      m+=1
    end
  end

  while(!eof(f) && (l=readline(f)) != "\$Elements")
  end

  l=readline(f)
  a=split(l, " ")
  blocks = parse(Int64, a[1])

  Triangles = Set(Array{Int64}[])
  Edges = Set(Array{Int64}[])
  for i=1:blocks
    l = readline(f)
    a = split(l, " ")
    elemDim = parse(Int64, a[1])
    elemEntitiy = parse(Int64, a[2])
    elemType = parse(Int64, a[3])
    elemsInBlock = parse(Int64, a[4])
    for j=1:elemsInBlock
      l = readline(f)
      a = split(l, " ")
      if(elemType == 1)
        push!(Edges, [tags[elemDim+1][elemEntitiy], parse(Int64, a[2]),
                      parse(Int64, a[3])])
      elseif(elemType == 2)
        push!(Triangles, [tags[elemDim+1][elemEntitiy], parse(Int64, a[2]),
                          parse(Int64, a[3]), parse(Int64, a[4])])
      else
        println("Not supported element tpye", elemType)
      end
    end
  end

  _Triangles = Array{Array{Int64,1},1}(undef,length(Triangles))
  _Edges = Array{Array{Int64,1},1}(undef,length(Edges))

  Volumes = Dict{Int64,Set{Int64}}()
  for (i,t) in enumerate(Triangles)
    _Triangles[i] = [NodeNumbering[n] for n in t[2:end]]
    if( !in(t[1],keys(Volumes)) )
      Volumes[t[1]] = Set{Int64}()
    end
    push!(Volumes[t[1]], i)
  end

  Boundaries = Dict{Int64, Boundary}()
  for (i,el) in enumerate(Edges)
    _Edges[i] = [NodeNumbering[n] for n in el[2:end]]
    if( !in(el[1],keys(Boundaries)) )
      Boundaries[el[1]] = Boundary(Set{Int64}(), Set{Int64}())
    end
    push!(Boundaries[el[1]].Edges, i)
    push!(Boundaries[el[1]].Nodes, NodeNumbering[el[2]])
    push!(Boundaries[el[1]].Nodes, NodeNumbering[el[3]])
  end

  return Mesh(nnodes, length(Triangles), length(Edges), Nodes, _Triangles, _Edges, Boundaries, Volumes)
end

function write_vtk_mesh(mesh::Mesh, file_name::String)
  points = zeros(Float64, length(mesh.Nodes[1]), length(mesh.Nodes))
  for (i,p) in enumerate(mesh.Nodes)
    points[:,i] = copy(p)
  end

  cells = Array{MeshCell,1}(undef, mesh.nelems)
  for (i,t) in enumerate(mesh.Triangles)
    cells[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, t)
  end

  return vtk_grid(file_name, points, cells)
end

function Jacobian(v1::Array{Float64,1}, v2::Array{Float64,1}, v3::Array{Float64,1})
  J = [v2-v1 v3-v1]
  return det(J), inv(J)'
end

function Jacobian(mesh::Mesh, elem::Int64)
  t = mesh.Triangles[elem]
  return Jacobian(mesh.Nodes[t[1]], mesh.Nodes[t[2]], mesh.Nodes[t[3]])
end

function EdgeJacobian(mesh::Mesh, elem::Int64)
  el = mesh.Edges[elem]
  return norm(mesh.Nodes[el[1]] - mesh.Nodes[el[2]])
end

function GetBoundaryNodes(mesh::Mesh, marker::Int64)
  return collect(mesh.Boundaries[marker].Nodes)
end

function getCellVolumes(mesh::Mesh)
  nelems = mesh.nelems
  v = zeros(nelems)
  for el=1:nelems
    nodes = mesh.Triangles[el]
    (detJ, J) = Jacobian(mesh, el)
    v[el] = sum(quadW)*detJ
  end
  return v
end
