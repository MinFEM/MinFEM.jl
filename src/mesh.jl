"""
$(TYPEDEF)

Abstract supertype for structs specifing regions of a domain, i.e. mesh.
Subtypes should contain at least Name, Nodes and Elements.
"""
abstract type Region end

"""
$(TYPEDEF)

Structure holding the name and sets of node and edge indices
for one particular physical boundary.

# Fields
$(TYPEDFIELDS)
"""
mutable struct Boundary <: Region
    "Unique physical name"
    Name::String
    
    "List of node indices"
    Nodes::Set{Int64}
    
    "List of element indices"
    Elements::Set{Int64}
end

Base.:(==)(a::Boundary, b::Boundary) = a.Nodes == b.Nodes &&
                                        a.Elements == b.Elements

"""
$(TYPEDEF)

Type holding the name and the set of element indices
for one particular physical domain.

# Fields
$(TYPEDFIELDS)
"""
mutable struct Domain <: Region
    "Unique physical name"
    Name::String

    "List of node indices"
    Nodes::Set{Int64}

    "List of element indices"
    Elements::Set{Int64}
end

Base.:(==)(a::Domain, b::Domain) = a.Elements == b.Elements

"""
$(TYPEDEF)

Type holding the associated physical tag and the set of elements
for one gmsh elementary entity.

# Fields
$(TYPEDFIELDS)
"""
mutable struct Entity
    "Unique physical id"
    PhysicalTag::Int64
    
    "List of element indices"
    Elements::Set{Int64}
end

Base.:(==)(a::Entity, b::Entity) = a.PhysicalTag == b.PhysicalTag && 
                                    a.Elements == b.Elements

"""
$(TYPEDEF)

Type for a triangular finite element mesh with volume and boundary markers.

# Fields
$(TYPEDFIELDS)
"""
mutable struct Mesh
    "Spatial dimension"
    d::Int64

    "Number of nodes"
    nnodes::Int64

    "Number of elements"
    nelems::Int64

    "Number of physical (marked) boundary elements"
    nboundelems::Int64

    "List of node coordinates"
    Nodes::Array{Array{Float64,1},1}

    "List of element node indices"
    Elements::Array{Array{Int64,1},1}

    "List of boundary element node indices"
    BoundaryElements::Array{Array{Int64,1},1}

    "List of parent elements to boundary elements"
    ParentElements::Array{Int64,1}

    "List of parent element local boundary to boundary elements"
    ParentBoundaries::Array{Int64,1}

    "Dictionary of marked boundaries"
    Boundaries::Dict{Int64,Boundary}

    "Dictionary of marked volume regions"
    Domains::Dict{Int64,Domain}

    "Dictionary of physical entities"
    Entities::Array{Dict{Int64,Entity},1}
end

Base.:(==)(a::Mesh, b::Mesh) = a.d == b.d && 
                                a.nnodes == b.nnodes &&
                                a.nelems == b.nelems && 
                                a.Nodes == b.Nodes &&
                                a.Elements == b.Elements && 
                                a.BoundaryElements == b.BoundaryElements && 
                                a.ParentElements == b.ParentElements && 
                                a.Boundaries == b.Boundaries && 
                                a.Domains == b.Domains &&
                                a.Entities == b.Entities

"""
$(TYPEDSIGNATURES)

Returns an n nodes quasi-uniform mesh for the 1D unit interval.
The left boundary node is denoted by 1001 and the right one by 1002.
"""
function unit_interval(n::Int64) :: Mesh
    if (n < 2)
        n = 2
    end
    d = 1
    nnodes = n
    nboundelems = 2
    nelems = n - 1
    h = 1 / (n - 1)

    Nodes = Array{Array{Float64,1},1}(undef, nnodes)
    BoundaryElements = Array{Array{Int64,1},1}(undef, nboundelems)
    Elements = Array{Array{Int64,1},1}(undef, nelems)
    Domains = Dict{Int64,Domain}()
    Boundaries = Dict{Int64,Boundary}()
    Boundaries[1001] = Boundary("left", Set{Int64}(), Set{Int64}())
    Boundaries[1002] = Boundary("right", Set{Int64}(), Set{Int64}())

    Nodes[1] = [0]
    Nodes[n] = [1]
    for k = 2:n-1
        Nodes[k] = [(k - 1) * h]
    end

    for i = 1:n-1
        Elements[i] = [i; i + 1]
    end

    push!(Boundaries[1001].Nodes, 1)
    push!(Boundaries[1002].Nodes, n)

    BoundaryElements[1] = [1]
    push!(Boundaries[1001].Elements, 1)

    BoundaryElements[2] = [n]
    push!(Boundaries[1002].Elements, 2)

    ParentElements::Array{Int64,1} = [1,nelems]
    ParentBoundaries::Array{Int64,1} = [1,2]


    Domains[10001] = Domain("domain", Set{Int64}(1:nnodes), Set{Int64}(1:nelems))

    Entities = [Dict{Int64,Entity}(), Dict{Int64,Entity}(),
                Dict{Int64,Entity}(), Dict{Int64,Entity}()]
    Entities[1][1] = Entity(1001, Set{Int64}(1))
    Entities[1][2] = Entity(1002, Set{Int64}(2))
    Entities[2][1] = Entity(10001, Set{Int64}(1:nelems))

    return Mesh(d, nnodes, nelems, nboundelems, Nodes, 
                Elements, BoundaryElements, ParentElements, ParentBoundaries,
                Boundaries, Domains, Entities)
end

"""
$(TYPEDSIGNATURES)

Returns a n-by-n quasi-uniform mesh for the 2D unit square.
The boundary indices are given in the order bottom, top, left, right from 1001 to 1004.
"""
function unit_square(n::Int64) :: Mesh
    if (n < 2)
        n = 2
    end
    d = 2
    nnodes = n * n
    nboundelems = 4 * (n - 1)
    nelems = 2 * (n - 1) * (n - 1)
    h = 1 / (n - 1)

    Nodes = Array{Array{Float64,1},1}(undef, nnodes)
    BoundaryElements = Array{Array{Int64,1},1}(undef, nboundelems)
    Elements = Array{Array{Int64,1},1}(undef, nelems)
    ParentElements = Array{Int64,1}(undef, nboundelems)
    ParentBoundaries = Array{Int64,1}(undef, nboundelems)
    Domains = Dict{Int64, Domain}()
    Boundaries = Dict{Int64,Boundary}()
    Boundaries[1001] = Boundary("bottom", Set{Int64}(), Set{Int64}())
    Boundaries[1002] = Boundary("top", Set{Int64}(), Set{Int64}())
    Boundaries[1003] = Boundary("left", Set{Int64}(), Set{Int64}())
    Boundaries[1004] = Boundary("right", Set{Int64}(), Set{Int64}())

    k = 1
    # bottom
    Nodes[k] = [0; 0]
    k+=1
    for j = 2:n-1
        Nodes[k] = [(j - 1) * h; 0]
        k += 1
    end
    Nodes[k] = [1; 0]
    k+=1

    # interior
    for i = 2:n-1
        Nodes[k] = [0; (i - 1) * h]
        k+=1
        for j = 2:n-1
            Nodes[k] = [(j - 1) * h; (i - 1) * h]
            k += 1
        end
        Nodes[k] = [1; (i - 1) * h]
        k+=1
    end

    # top
    Nodes[k] = [0; 1]
    k+=1
    for j = 2:n-1
        Nodes[k] = [(j - 1) * h; 1]
        k += 1
    end
    Nodes[k] = [1; 1]

    k = 1
    for i = 1:n-1
        for j = 1:n-1
            Elements[k] = [(i - 1) * n + j; (i - 1) * n + j + 1; i * n + j]
            k = k + 1

            Elements[k] = [(i - 1) * n + j + 1; i * n + j + 1; i * n + j]
            k = k + 1
        end
    end

    for j = 1:n
        push!(Boundaries[1001].Nodes, j)
        push!(Boundaries[1002].Nodes, j + (n - 1) * n)
        push!(Boundaries[1003].Nodes, (j - 1) * n + 1)
        push!(Boundaries[1004].Nodes, (j - 1) * n + n)
    end

    k = 1
    for j = 1:n-1
        BoundaryElements[k] = [j, j + 1]
        ParentElements[k] = j * 2 - 1
        ParentBoundaries[k] = 1
        push!(Boundaries[1001].Elements, k)
        k = k + 1
    end

    for j = 1:n-1
        BoundaryElements[k] = [j + 1 + (n - 1) * n, j + (n - 1) * n]
        ParentElements[k] = j * 2 + 2 * (n-1) * (n-2)
        ParentBoundaries[k] = 3
        push!(Boundaries[1002].Elements, k)
        k = k + 1
    end

    for j = 1:n-1
        BoundaryElements[k] = [(j - 1) * n + 1, j * n + 1]
        ParentElements[k] = (j-1) * 2 * (n-1) + 1
        ParentBoundaries[k] = 2
        push!(Boundaries[1003].Elements, k)
        k = k + 1
    end
    
    for j = 1:n-1
        BoundaryElements[k] = [(j - 1) * n + n, j * n + n]
        ParentElements[k] = (j-1) * 2 * (n-1) + 2 * (n-1)
        ParentBoundaries[k] = 1
        push!(Boundaries[1004].Elements, k)
        k = k + 1
    end

    Domains[10001] = Domain("domain", Set{Int64}(1:nnodes), Set{Int64}(1:nelems))

    Entities = [Dict{Int64,Entity}(), Dict{Int64,Entity}(), 
                Dict{Int64,Entity}(), Dict{Int64,Entity}()]
    Entities[2][1] = Entity(1001, copy(Boundaries[1001].Elements))
    Entities[2][2] = Entity(1002, copy(Boundaries[1002].Elements))
    Entities[2][3] = Entity(1003, copy(Boundaries[1003].Elements))
    Entities[2][4] = Entity(1004, copy(Boundaries[1004].Elements))
    Entities[3][1] = Entity(10001, Set{Int64}(1:nelems))

    return Mesh(d, nnodes, nelems, nboundelems, Nodes, 
                Elements, BoundaryElements, ParentElements, ParentBoundaries,
                Boundaries, Domains, Entities)
end


"""
$(TYPEDSIGNATURES)

Returns a mesh imported from a gmsh file of version v1, v2 or v4.
"""
function import_mesh(fileName::String) :: Mesh
    f = open(fileName)

    v1::Bool = false
    while !eof(f)
        l=readline(f)
        if l == "\$NOD"
            v1 = true
            break
        elseif l == "\$MeshFormat"
            break
        end
    end

    if eof(f)
        throw(ErrorException("Mesh file corrupted."))
    elseif v1
        mesh = import_mesh1(f)
    else
        l=readline(f)
        a = split(l, " ")
        version = parse(Float64, a[1])
        if version >= 2.0 && version < 3.0
            mesh = import_mesh2(f)
        elseif version >= 4.1 && version < 5.0
            mesh = import_mesh4(f)
        else
            throw(ErrorException("Unsupported mesh format: " *
                                    "$(string(parse(Float64, a[1]))).\n" *
                                    "msh4 is recommended."))
        end
    end

    close(f)
    return mesh
end

"""
$(TYPEDSIGNATURES)

Returns mesh by continuing import process started by `import_mesh()`
for gmsh files of version 1. 
"""
function import_mesh1(f::IOStream)
    l=readline(f)
    nnodes = parse(Int64, l)
    _Nodes = Array{Array{Float64, 1}, 1}(undef, nnodes)
    for i=1:nnodes
        l = readline(f)
        a = split(l, " ")
        _Nodes[i] = parse.(Float64, a[2:4])
    end

    while(!eof(f) && (l=readline(f)) != "\$ELM")
    end
    l=readline(f)
    nelems = parse(Int64, l)
  
    _Elements = Dict{Int64, Array{Array{Int64,1},1}}()
    for i = 1:nelems
        l = readline(f)
        a = split(l, " ")

        elemDim = gmsh_dimfromtype(parse(Int64, a[2]))
        elemTag = parse(Int64, a[3])
        elemEntitity = parse(Int64, a[4])
        elemNodes = parse(Int64, a[5])

        if !in(elemDim, keys(_Elements))
            _Elements[elemDim] = Array{Array{Int64,1},1}(undef,0)
        end

        val = [2, elemTag, elemEntitity]
        append!(val, parse.(Int64, a[6:(5+elemNodes)]))
        append!(_Elements[elemDim], [val])
    end

    d::Int64 = maximum(collect(keys(_Elements)))

    Nodes = Array{Array{Float64, 1}, 1}(undef, nnodes)
    for i = 1:nnodes
        Nodes[i] = copy(_Nodes[i][1:d])
    end
  
    Elements = Array{Array{Int64,1},1}(undef,length(_Elements[d]))
    BoundaryElements = Array{Array{Int64,1},1}(undef,length(_Elements[d-1]))
    ParentElements = Array{Int64,1}(undef,length(_Elements[d-1]))
    ParentBoundaries = Array{Int64,1}(undef,length(_Elements[d-1]))
    Domains = Dict{Int64, Domain}()
    Boundaries = Dict{Int64, Boundary}()
    Entities = [Dict{Int64,Entity}(), Dict{Int64,Entity}(), 
                Dict{Int64,Entity}(), Dict{Int64,Entity}()]

    for (i,el) in enumerate(_Elements[d])
        Elements[i] = copy(el[4:end])
        if !in(el[2], keys(Domains))
            Domains[el[2]] = Domain("", Set{Int64}(), Set{Int64}())
        end
        push!(Domains[el[2]].Elements, i)
        for node in el[4:end]
            push!(Domains[el[2]].Nodes, node)
        end

        if !in(el[3],keys(Entities[d+1]))
            Entities[d+1][el[3]] = Entity(el[2], Set{Int64}())
        end
        push!(Entities[d+1][el[3]].Elements, i)
    end
  
    for (i,el) in enumerate(_Elements[d-1])
        _boundarynodes = copy(el[4:end])

        ParentElements[i] = findfirst(x -> issubset(_boundarynodes, x), Elements)
        ParentBoundaries[i] = parentboundary(
            _boundarynodes, 
            Elements[ParentElements[i]]
        )
        BoundaryElements[i] = sort_boundaryelement(
            _boundarynodes,
            Elements[ParentElements[i]]
        )
 
        if !in(el[2],keys(Boundaries))
            Boundaries[el[2]] = Boundary("", Set{Int64}(), Set{Int64}())
        end
        push!(Boundaries[el[2]].Elements, i)
        for node in _boundarynodes
            push!(Boundaries[el[2]].Nodes, node)
        end

        if !in(el[3],keys(Entities[d]))
            Entities[d][el[3]] = Entity(el[2], Set{Int64}())
        end
        push!(Entities[d][el[3]].Elements, i)
    end
  
    return Mesh(
        d,
        nnodes,
        length(Elements),
        length(BoundaryElements),
        Nodes, 
        Elements,
        BoundaryElements,
        ParentElements,
        ParentBoundaries,
        Boundaries,
        Domains,
        Entities
    )
  end

  """
$(TYPEDSIGNATURES)

Returns mesh by continuing import process started by `import_mesh()`
for gmsh files of version 2. 
"""
function import_mesh2(f::IOStream)
    physicalNames = [Dict{Int64,String}(), Dict{Int64,String}(),
                        Dict{Int64,String}(), Dict{Int64,String}()]

    while !eof(f)
        l=readline(f)
        if l == "\$PhysicalNames"
            l = readline(f)
            a = split(l, " ")
            numPhysicalNames = parse(Int64, a[1])
            
            for r = 1:numPhysicalNames
                l = readline(f)
                a = split(l, " ")
                id1 = parse(Int64, a[1]) + 1
                id2 = parse(Int64, a[2])
                physicalNames[id1][id2] = chop(a[3], head = 1, tail = 1)
            end
            # No break and run until other case

        elseif l == "\$Nodes"
            break
        end
    end

    l=readline(f)
    nnodes = parse(Int64, l)
    _Nodes = Array{Array{Float64, 1}, 1}(undef, nnodes)
    for k = 1:nnodes
        l = readline(f)
        a = split(l, " ")
        _Nodes[k] = parse.(Float64, a[2:4])
    end

    while(!eof(f) && (l=readline(f)) != "\$Elements")
    end
    l=readline(f)
    nelems = parse(Int64, l)
  
    _Elements = [Array{Array{Int64,1},1}(), Array{Array{Int64,1},1}(),
                    Array{Array{Int64,1},1}(), Array{Array{Int64,1},1}()]
    for i = 1:nelems
        l = readline(f)
        a = split(l, " ")

        elemDim = gmsh_dimfromtype(parse(Int64, a[2]))
        val = parse.(Int64, a[3:(6+elemDim)])
        append!(_Elements[elemDim+1], [val])
    end

    # Find dimension
    d::Int64 = 0
    for j = 2:4
        if !isempty(_Elements[j])
            d = j-1
        end
    end
    if d == 0
        throw(ErrorException("File does not contain a meshed domain."))
    end

    Nodes = Array{Array{Float64, 1}, 1}(undef, nnodes)
    for k = 1:nnodes
        Nodes[k] = _Nodes[k][1:d]
    end
  
    Elements = Array{Array{Int64,1},1}(undef, length(_Elements[d+1]))
    BoundaryElements = Array{Array{Int64,1},1}(undef, length(_Elements[d]))
    ParentElements = Array{Int64,1}(undef, length(BoundaryElements))
    ParentBoundaries = Array{Int64,1}(undef, length(BoundaryElements))
    Domains = Dict{Int64, Domain}()
    Boundaries = Dict{Int64, Boundary}()
    Entities = [Dict{Int64,Entity}(), Dict{Int64,Entity}(),
                Dict{Int64,Entity}(), Dict{Int64,Entity}()]

    for (i,el) in enumerate(_Elements[d+1])
        Elements[i] = copy(el[4:end])
        if !in(el[2], keys(Domains))
            Domains[el[2]] = Domain(
                get!(physicalNames[d+1],el[2],""),
                Set{Int64}(),
                Set{Int64}()
            )
        end
        push!(Domains[el[2]].Elements, i)
        for node in el[4:end]
            push!(Domains[el[2]].Nodes, node)
        end

        if !in(el[3],keys(Entities[d+1]))
            Entities[d+1][el[3]] = Entity(el[2], Set{Int64}())
        end
        push!(Entities[d+1][el[3]].Elements, i)
    end
  
    for (i,el) in enumerate(_Elements[d])
        _boundarynodes = copy(el[4:end])

        ParentElements[i] = findfirst(x -> issubset(_boundarynodes, x), Elements)
        ParentBoundaries[i] = parentboundary(
            _boundarynodes, 
            Elements[ParentElements[i]]
        )
        BoundaryElements[i] = sort_boundaryelement(
            _boundarynodes,
            Elements[ParentElements[i]]
        )

        if !in(el[2],keys(Boundaries))
            Boundaries[el[2]] = Boundary(get!(physicalNames[d],el[2],""),
                                            Set{Int64}(), Set{Int64}())
        end
        push!(Boundaries[el[2]].Elements, i)
        for node in _boundarynodes
            push!(Boundaries[el[2]].Nodes, node)
        end

        if !in(el[3],keys(Entities[d]))
            Entities[d][el[3]] = Entity(el[2], Set{Int64}())
        end
        push!(Entities[d][el[3]].Elements, i)
    end
  
    return Mesh(
        d,
        nnodes,
        length(Elements),
        length(BoundaryElements),
        Nodes, 
        Elements,
        BoundaryElements,
        ParentElements,
        ParentBoundaries,
        Boundaries,
        Domains,
        Entities
    )
end

"""
$(TYPEDSIGNATURES)

Returns mesh by continuing import process started by `import_mesh()` 
for gmsh files of version 4. 
"""
function import_mesh4(f::IOStream)
    physicalNames = [Dict{Int64,String}(), Dict{Int64,String}(),
                        Dict{Int64,String}(), Dict{Int64,String}()]

    while !eof(f)
        l=readline(f)
        if l == "\$PhysicalNames"
            l = readline(f)
            a = split(l, " ")
            numPhysicalNames = parse(Int64, a[1])
            
            for r = 1:numPhysicalNames
                l = readline(f)
                a = split(l, " ")
                id1 = parse(Int64, a[1]) + 1
                id2 = parse(Int64, a[2])
                physicalNames[id1][id2] = chop(a[3], head = 1, tail = 1)
            end
            # No break and run until other case

        elseif l == "\$Entities"
            break
        end
    end

    # 0,1,2,3 dimensional entity tags
    Entities = [Dict{Int64,Entity}(), Dict{Int64,Entity}(),
                Dict{Int64,Entity}(), Dict{Int64,Entity}()]
    l = readline(f)
    a = split(l, " ")
    numTags = [parse(Int64, a[1]), parse(Int64, a[2]), 
                parse(Int64, a[3]), parse(Int64, a[4])]

    # Point tags
    for r = 1:numTags[1]
        l = readline(f)
        a = split(l, " ")
        if (parse(Int64, a[5]) != 0) # Check if physical or not and parse name
            Entities[1][parse(Int64, a[1])] = Entity(parse(Int64, a[6]), Set{Int64}())
        end
    end

    # Curve, Surface, Volume tags
    d::Int64 = 0
    for j = 2:4
        for r = 1:numTags[j]
            l = readline(f)
            a = split(l, " ")
            if (parse(Int64, a[8]) != 0) # Check if physical or not
                Entities[j][parse(Int64, a[1])] = Entity(parse(Int64, a[9]), Set{Int64}())
            end
        end
        if !isempty(Entities[j])
            d = j-1
        end
    end

    while (!eof(f) && (l = readline(f)) != "\$Nodes")
    end
    if eof(f)
        throw(ErrorException("File does not contain a meshed domain."))
    end

    l = readline(f)
    a = split(l, " ")
    blocks = parse(Int64, a[1])
    nnodes = parse(Int64, a[2])

    Nodes = Array{Array{Float64,1},1}(undef, nnodes)
    NodeNumbering = Dict{Int64,Int64}()
    n = 1
    m = 1
    for i = 1:blocks
        l = readline(f)
        a = split(l, " ")
        nodesInBlock = parse(Int64, a[4])

        # The node numbers
        for j = 1:nodesInBlock
            l = readline(f)
            a = split(l, " ")
            NodeNumbering[parse(Int64, a[1])] = n
            n += 1
        end

        # The actual coordinates
        for j = 1:nodesInBlock
            l = readline(f)
            a = split(l, " ")
            Nodes[m] = parse.(Float64, a[1:d])
            m += 1
        end
    end

    while (!eof(f) && (l = readline(f)) != "\$Elements")
    end

    l = readline(f)
    a = split(l, " ")
    blocks = parse(Int64, a[1])

    _Elements = Array{Array{Int64,1},1}()
    _BoundaryElements = Array{Array{Int64,1},1}()
    boundaryElementType = gmsh_typefromdim(d-1)
    elementType = gmsh_typefromdim(d)
    for i = 1:blocks
        l = readline(f)
        a = split(l, " ")
        elemDim = parse(Int64, a[1])
        elemEntitity = parse(Int64, a[2])
        elemType = parse(Int64, a[3])
        elemsInBlock = parse(Int64, a[4])
        elemTag = Entities[elemDim+1][elemEntitity].PhysicalTag
        
        for j = 1:elemsInBlock
            l = readline(f)
            a = split(l, " ")
            val = [2, elemTag, elemEntitity]
            append!(val, parse.(Int64, a[2:(2+elemDim)]))

            if elemType == boundaryElementType
                append!(_BoundaryElements, [val])
            elseif elementType == elementType
                append!(_Elements, [val])
            else
                println("Not supported element tpye $elemType for $d-dimensional mesh.")
            end
        end
    end

    Elements = Array{Array{Int64,1},1}(undef, length(_Elements))
    BoundaryElements = Array{Array{Int64,1},1}(undef, length(_BoundaryElements))
    ParentElements = Array{Int64,1}(undef,length(_BoundaryElements))
    ParentBoundaries = Array{Int64,1}(undef,length(_BoundaryElements))
    Domains = Dict{Int64,Domain}()
    Boundaries = Dict{Int64,Boundary}()

    for (i, el) in enumerate(_Elements)
        Elements[i] = [NodeNumbering[n] for n in el[4:end]]
        if (!in(el[2], keys(Domains)))
            Domains[el[2]] = Domain(
                get!(physicalNames[d+1],el[2],""),
                Set{Int64}(),
                Set{Int64}()
            )
        end
        push!(Domains[el[2]].Elements, i)
        for node in el[4:end]
            push!(Domains[el[2]].Nodes, NodeNumbering[node])
        end
        
        push!(Entities[d+1][el[3]].Elements, i)
    end

    for (i, el) in enumerate(_BoundaryElements)
        _boundarynodes = [NodeNumbering[n] for n in el[4:end]]

        ParentElements[i] = findfirst(x -> issubset(_boundarynodes, x), Elements)
        ParentBoundaries[i] = parentboundary(
            _boundarynodes, 
            Elements[ParentElements[i]]
        )
        BoundaryElements[i] = sort_boundaryelement(
            _boundarynodes,
            Elements[ParentElements[i]]
        )

        if (!in(el[2], keys(Boundaries)))
            Boundaries[el[2]] = Boundary(get!(physicalNames[d],el[2],""),
                                            Set{Int64}(), Set{Int64}())
        end
        
        push!(Boundaries[el[2]].Elements, i)
        for n in _boundarynodes
            push!(Boundaries[el[2]].Nodes, n)
        end

        push!(Entities[d][el[3]].Elements, i)
    end

    return Mesh(
        d,
        nnodes,
        length(Elements),
        length(BoundaryElements),
        Nodes, 
        Elements,
        BoundaryElements,
        ParentElements,
        ParentBoundaries,
        Boundaries,
        Domains,
        Entities
    )
end

"""
$(TYPEDSIGNATURES)

Exports a mesh to a gmsh file of version v2.
"""
function export_mesh(mesh::Mesh, fileName::String)
    if !occursin(".msh",fileName)
        fileName = fileName * ".msh2"
    end
    f = open(fileName, "w")
    
    write(f, "\$MeshFormat\n")
    write(f, "2.2 0 8\n")
    write(f, "\$EndMeshFormat\n")

    nPhysicalNames = length(findall(x -> x.Name != "", mesh.Boundaries)) +
                        length(findall(x -> x.Name != "", mesh.Domains))
    if nPhysicalNames > 0
        write(f, "\$PhysicalNames\n")
        write(f, "$(length(mesh.Boundaries)+length(mesh.Domains))\n")
        for (key,val) in sort(collect(pairs(mesh.Boundaries)), by=x->x[1])
            if val.Name != ""
                write(f, "$(mesh.d-1) $key \"$(val.Name)\"\n")
            end
        end
        for (key,val) in sort(collect(pairs(mesh.Domains)), by=x->x[1])
            if val.Name != ""
                write(f, "$(mesh.d) $key \"$(val.Name)\"\n")
            end
        end
        write(f, "\$EndPhysicalNames\n")
    end
    
    write(f, "\$Nodes\n")
    write(f, "$(mesh.nnodes)\n")
    for i = 1:mesh.nnodes
        write(f, "$i")
        for j = 1:mesh.d
            write(f," $(mesh.Nodes[i][j])")
        end
        for j = mesh.d:2
            write(f," 0")
        end
        write(f," \n")
    end
    write(f, "\$EndNodes\n")
    
    write(f, "\$Elements\n")
    write(f, "$(mesh.nelems+mesh.nboundelems)\n")
    k = 1
    for i in [mesh.d, mesh.d+1]
        for (key,val) in sort(collect(pairs(mesh.Entities[i])), by=x->x[1])
            for el in sort(collect(val.Elements))
                if i == mesh.d
                    elem = mesh.BoundaryElements[el]
                else
                    elem = mesh.Elements[el]
                end
                gmshElemType = gmsh_typefromdim(i-1) 
                write(f, "$k $gmshElemType 2 $(val.PhysicalTag) $key")
                for j in eachindex(elem)
                    write(f, " $(elem[j])")
                end
                write(f,"\n")
                k += 1
            end
        end
    end
    write(f, "\$EndElements\n")

    close(f)
end

"""
$(TYPEDSIGNATURES)
    
Returns gmsh elementary element type for given dimension.
"""
function gmsh_typefromdim(d::Int64)
    if d == 0
        return 15 # 1-node point
    elseif d == 1
        return 1 # 2-node line
    elseif d == 2
        return 2 # 3-node triangle
    elseif d == 3
        return 4 # 4-node tetrahedron
    else
        throw(DomainError(d, "Unsupported dimension. " *
                                "Only dimensions 0,1,2 and 3 are supported."))
    end
end

"""
$(TYPEDSIGNATURES)
    
Returns dimension of a gmsh elementary entitiy type.
"""
function gmsh_dimfromtype(t::Int64)
    if t == 15
        return 0 # 1-node point
    elseif t == 1
        return 1 # 2-node line
    elseif t == 2
        return 2 # 3-node triangle
    elseif t == 4
        return 3 # 4-node tetrahedron
    else
        throw(DomainError(t, "Unsupported element type. " *
                                "Only first order tetrahedral and corresponing" *
                                "lower dimensional types are supported."))
    end
end

"""
$(TYPEDSIGNATURES)
    
Returns index of the boundary of the parent element spanned by the nodes.
"""
function parentboundary(nodes::Vector{Int64}, parentNodes::Vector{Int64})
    dim = length(nodes)

    localindex = zeros(Int64, dim)
    for j = 1:dim
        localindex[j] = findfirst(x -> x == nodes[j], parentNodes)
    end

    if dim == 1
        return findfirst(x -> x == nodes[1], parentNodes)
    elseif dim == 2
        if issetequal(localindex, [1,2])
            return 1
        elseif issetequal(localindex, [1,3])
            return 2
        elseif issetequal(localindex, [2,3])
            return 3
        else
            throw(DomainError(nodes, "Nodes are not subset of parent."))
        end
    elseif dim == 3
        if issetequal(localindex, [1,2,3])
            return 1
        elseif issetequal(localindex, [1,2,4])
            return 2
        elseif issetequal(localindex, [1,3,4])
            return 3
        elseif issetequal(localindex, [2,3,4])
            return 4
        else
            throw(DomainError(nodes, "Nodes are not subset of parent."))
        end
    end
end

"""
$(TYPEDSIGNATURES)
    
Sort nodes in boundary element corresponding to the order they occur in the parent element.
Can be used to determine orientation of the edge.
"""
function sort_boundaryelement(nodes::Array{Int64}, parent::Array{Int64})
    sorted = Array{Int64,1}(undef,length(nodes))
    k = 1
    for node in parent
        if node in nodes
            sorted[k] = node
            k = k + 1
        end
    end

    return sorted
end

"""
$(TYPEDSIGNATURES)
    
Updates given mesh by shifting all nodes to new coordinates c.
"""
function update_mesh!(mesh::Mesh, c::Array{Array{Float64,1},1})
    if length(c) != mesh.nnodes
        throw(ArgumentError("Deformation vector does not have matching length."))
    end

    mesh.Nodes = deepcopy(c)
end

"""
    deform_mesh!(mesh::Mesh, v::AbstractVector{Float64}; t::Float64=1.0)
    
Deforms given mesh by shifting all nodes according to the vector field v 
scaled by the stepsize t.
"""
function deform_mesh!(mesh::Mesh, v::AbstractVector{Float64}; t::Float64=1.0)
    if length(v) != mesh.nnodes * mesh.d
        throw(ArgumentError("Deformation vector does not have matching length."))
    end

    for i = 1:mesh.nnodes
        mesh.Nodes[i] += t * v[mesh.d*(i-1)+1 : mesh.d*i]
    end
end

"""
    deform_mesh(mesh::Mesh, v::AbstractVector{Float64}; t::Float64=1.0) -> Mesh
    
Returns copy of mesh deformed by `deform_mesh!`.
"""
function deform_mesh(mesh::Mesh, v::AbstractVector{Float64}; t::Float64=1.0)
    _mesh = deepcopy(mesh)
    deform_mesh!(_mesh, v, t=t)
    return _mesh
end

"""
$(TYPEDSIGNATURES)

Returns set of all or specified physical boundaries of the mesh.
"""
function select_boundaries(mesh::Mesh, args...)
    Boundaries = Array{Boundary,1}()
    
    for id in args
        if in(id, keys(mesh.Boundaries))
            push!(Boundaries, mesh.Boundaries[id])
        else
            throw(ArgumentError("Physical boundary id $(id) not valid integer."))
        end
    end

    if length(args) == 0
        for (_, bnd) in mesh.Boundaries
            push!(Boundaries, bnd)
        end
    end

    return Set(Boundaries)
end

"""
 $(TYPEDSIGNATURES)

Returns set of all or specified physical boundaries of the mesh.
"""
function select_domains(mesh::Mesh, args...)
    Domains = Array{Domain,1}()
    
    for id in args
        if in(id, keys(mesh.Domains))
            push!(Domains, mesh.Domains[id])
        else
            throw(ArgumentError("Physical domain id $(id) not valid integer."))
        end
    end

    if length(args) == 0
        for (_, domain) in mesh.Domains
            push!(Domains, domain)
        end
    end

    return Set(Domains)
end

"""
$(TYPEDSIGNATURES)

Returns set of node ids in set of physical boundaries.
"""
function extract_nodes(boundaries::Set{Boundary})
    nodes = Set{Int64}()
    for boundary in boundaries
        union!(nodes, boundary.Nodes)
    end
    return nodes
end

"""
$(TYPEDSIGNATURES)

Returns set of node ids in set of physical domains.
"""
function extract_nodes(domains::Set{Domain})
    nodes = Set{Int64}()
    for domain in domains
        union!(nodes, domain.Nodes)
    end
    return nodes
end

"""
$(TYPEDSIGNATURES)

Returns set of boundary element ids in set of physical boundaries.
"""
function extract_elements(boundaries::Set{Boundary})
    elements = Set{Int64}()
    for boundary in boundaries
        union!(elements, boundary.Elements)
    end
    return elements
end

"""
$(TYPEDSIGNATURES)

Returns set of boundary element ids in set of physical domains.
"""
function extract_elements(domains::Set{Domain})
    elements = Set{Int64}()
    for domain in domains
        union!(elements, domain.Elements)
    end
    return elements
end

"""
    evaluate_mesh_function(
        mesh::Mesh,
        f::Function; 
        region::Set{Int64} = Set{Int64}(),
        qdim::Int64 = 1
    ) -> Vector{Float64}

Returns evaluation of a given function object f on all or specified nodes of the mesh.
Can be either called with set of physical boundaries or directly with a set of nodes 
when given with keyword argument region.
"""
function evaluate_mesh_function(
    mesh::Mesh,
    f::Function; 
    region::Set{Int64} = Set{Int64}(),
    qdim::Int64 = 1
)
    v = zeros(Float64, qdim * mesh.nnodes)
    nodes = 1:mesh.nnodes

    if !isempty(region)
        nodes = collect(region)
    end

    if qdim == 1
        for i in nodes
            v[i] = f(mesh.Nodes[i])
        end
    else
        for i in nodes
            v[qdim*(i-1)+1 : qdim*i] = f(mesh.Nodes[i])
        end
    end

    return v
end

"""
    evaluate_mesh_function(
        mesh::Mesh,
        f::Function,
        region::Set{Boundary};
        qdim::Int64 = 1
    ) -> Vector{Float64}

Same as previous $(FUNCTIONNAME)(...) but takes `Set{Boundary}` as mandatory argument
for the relevant region, which replaces the keyword argument using `Set{Int64}`.
Extracts nodes from the boundary and then passes them to the base function. 
"""
function evaluate_mesh_function(
    mesh::Mesh,
    f::Function,
    region::Set{Boundary};
    qdim::Int64 = 1
)
    evaluate_mesh_function(mesh, f, region = extract_nodes(region), qdim = qdim)
end

"""
    evaluate_function(
        f::Vector{Float64},
        node::Int64;
        qdim::Int64 = 1
    ) -> Vector{Float64}

Evaluates discrete function give as FEM coefficient vector at a given node.
Resulting vector has the lenght of the number of components of f. 
"""
function evaluate_function(
    f::Vector{Float64},
    node::Int64;
    qdim::Int64 = 1
)
    return f[qdim*(node-1)+1 : qdim*node]
end

"""
    evaluate_function(
        f::Vector{Float64},
        region::Set{Int64};
        qdim::Int64 = 1
    ) -> Vector{Float64}

Evaluates discrete function give as FEM coefficient vector at all nodes in a given region.
Resulting vector has the lenght of the full coefficient vector with the entries
corresponding to non-evaluated nodes being zero. 
"""
function evaluate_function(
    f::Vector{Float64},
    region::Set{Int64};
    qdim::Int64 = 1
)
    l = length(f)
    v = zeros(Float64, l)

    for i in region
        v[qdim*(i-1)+1 : qdim*i] = f[qdim*(i-1)+1 : qdim*i]
    end

    return v
end

"""
$(TYPEDSIGNATURES)
    
Returns transformation matrix (jacobian) of the mapping 
from the FEM reference element to an element spanned by nodes at the given coordinates.

Commonly the coordinates correspond to an element in a mesh, but not necessarily have to.
"""
function base_jacobian(coords::Vector{Vector{Float64}})
    d = length(coords) - 1
    
    J = zeros(Float64, d, d)
    for i = 1:d
        J[:,i] = coords[i+1] - coords[1] 
    end

    return J
end

"""
$(TYPEDSIGNATURES)
    
Same as previous `$(FUNCTIONNAME)(...)`, but takes a mesh and set of nodes as arguments.
Extracts coordinates from the mesh and passes them to the base function.

Commonly the coordinates correspond to one element in a mesh, but not necessarily have to.
"""
function base_jacobian(mesh::Mesh, nodes::Array{Int64,1})
    return base_jacobian(mesh.Nodes[nodes])
end

"""
$(TYPEDSIGNATURES)
    
Same as previous `$(FUNCTIONNAME)(...)`, but takes a mesh and an element id as arguments.
Extracts coordinates of the support nodes from the mesh
and passes them to the base function.
"""
function base_jacobian(mesh::Mesh, element::Int64)
    return base_jacobian(mesh.Nodes[mesh.Elements[element]])
end

"""
$(TYPEDSIGNATURES)
    
Returns determinant (i.e. element weight) and inverse transposed of the jacobian 
of an FEM element spanned by nodes at the given coordinates.

Commonly the coordinates correspond to one element in a mesh, but not necessarily have to.
"""
function jacobian(
    coords::Vector{Vector{Float64}}
) :: Tuple{Float64, AbstractMatrix{Float64}}
    J = base_jacobian(coords)
    return det(J), inv(J)'
end

"""
$(TYPEDSIGNATURES)
    
Same as previous `$(FUNCTIONNAME)(...)`, but takes a mesh and set of nodes as arguments.
Extracts coordinates from the mesh and passes them to the base function.

Commonly the coordinates correspond to one element in a mesh, but not necessarily have to.
"""
function jacobian(mesh::Mesh, nodes::Array{Int64,1})
    return jacobian(mesh.Nodes[nodes])
end

"""
$(TYPEDSIGNATURES)
    
Same as previous `$(FUNCTIONNAME)(...)`, but takes a mesh and an element id as arguments.
Extracts coordinates of the support nodes from the mesh
and passes them to the base function.
"""
function jacobian(mesh::Mesh, element::Int64)
    return jacobian(mesh, mesh.Elements[element])
end

"""
$(TYPEDSIGNATURES)
    
Returns determinant of the jacobian (i.e. element weight) of an FEM boundary element
(in d-1 dimensions) spanned spanned by nodes at the given coordinates.

Commonly the coordinates correspond to an element in a mesh, but not necessarily have to.
"""
function jacobian_boundary(coords::Vector{Vector{Float64}})
    d = length(coords)
    
    if d == 3
        return norm(cross(coords[3] - coords[1], coords[2] - coords[1]), 2)
    elseif d == 2
        return norm(coords[2] - coords[1], 2)
    elseif d == 1
        return length(coords)
    end
end

"""
$(TYPEDSIGNATURES)
    
Same as previous `$(FUNCTIONNAME)(...)`, but takes a mesh and set of nodes as arguments.
Extracts coordinates from the mesh and passes them to the base function.

Commonly the coordinates correspond to one element in a mesh, but not necessarily have to.
"""
function jacobian_boundary(mesh::Mesh, nodes::Array{Int64,1})
    return jacobian_boundary(mesh.Nodes[nodes])
end

"""
$(TYPEDSIGNATURES)
    
Same as previous `$(FUNCTIONNAME)(...)`, but takes a mesh and an element id as arguments.
Extracts coordinates of the support nodes from the mesh
and passes them to the base function.
"""
function jacobian_boundary(mesh::Mesh, element::Int64)
    return jacobian_boundary(mesh.Nodes[mesh.BoundaryElements[element]])
end

"""
$(TYPEDSIGNATURES)
    
Returns the outer normal vector at the ii-th boundary of the
dim-dimensional reference element. 
"""
function outernormalvector(
    dim::Int64,
    ii::Int64
)
    if ii == dim+1
        return ones(Float64, dim) ./ sqrt(dim)
    else
        eta = zeros(Float64, dim)
        eta[dim+1-ii] = -1
        return eta
    end
end

"""
$(TYPEDSIGNATURES)
    
Returns the outer normal vector at a boundary element of the given mesh using a
pre-computed jacobian transformation matrix of the parent element.
"""
function outernormalvector(
    mesh::Mesh,
    boundaryElement::Int64,
    J::AbstractMatrix{Float64}
)
    refNormal = outernormalvector(mesh.d, mesh.ParentBoundaries[boundaryElement])
    
    mesh.d == 1 && return refNormal

    orth = J * refNormal
    return orth ./ norm(orth,2)
end

"""
$(TYPEDSIGNATURES)
    
Returns the outer normal vector at the given boundary element of the given mesh.
"""
function outernormalvector(
    mesh::Mesh,
    boundaryElement::Int64
)
    _, J = jacobian(mesh, mesh.Elements[mesh.ParentElements[boundaryElement]])
    return outernormalvector(mesh, boundaryElement, J)
end

"""
    outernormalvector(
        mesh::Mesh;
        boundaryElements::Set{Int64} = Set{Int64}()
    ) -> Vector{Float64}
    
Returns coefficient vector of outer normal vectors at all or specified boundary elements
of the given mesh.
"""
function outernormalvector(
    mesh::Mesh;
    boundaryElements::Set{Int64} = Set{Int64}()
)
    if isempty(boundaryElements)
        boundaryElements = Set{Int64}(1 : mesh.nboundelems)
    end
    
    eta = zeros(Float64, mesh.d * mesh.nboundelems)
    
    for bel in boundaryElements
        eta[mesh.d*(bel-1)+1 : mesh.d*bel] += outernormalvector(mesh, bel)
    end

    return eta
end

"""
$(TYPEDSIGNATURES)

Returns volume of the d-dimensional reference element. 
"""
function elementvolume(d::Int64)
    return 1 / factorial(d)
end

"""
$(TYPEDSIGNATURES)

Returns volume of the given element in the given mesh. 
"""
function elementvolume(mesh::Mesh, element::Int64)
    detJ = det(base_jacobian(mesh, element))
    return detJ * elementvolume(mesh.d)
end

"""
$(TYPEDSIGNATURES)

Returns vector of volumes of all elements in the given mesh.
"""
function elementvolume(mesh::Mesh)
    v = zeros(mesh.nelems)
    ref_vol = elementvolume(mesh.d)

    for el in eachindex(v)
        detJ = det(base_jacobian(mesh, el))
        v[el] = detJ * ref_vol
    end

    return v
end

"""
$(TYPEDSIGNATURES)

Returns volume of the (d-1)-dimensional boundary reference element. 
"""
function elementvolume_boundary(d::Int64)
    return elementvolume(d-1)
end

"""
$(TYPEDSIGNATURES)

Returns volume of the given boundary element in the given mesh. 
"""
function elementvolume_boundary(mesh::Mesh, element::Int64)
    detJ = jacobian_boundary(mesh, element)
    return detJ * elementvolume(mesh.d-1)
end

"""
$(TYPEDSIGNATURES)

Returns vector of volumes of all boundary elements in the given mesh.
"""
function elementvolume_boundary(mesh::Mesh)
    v = zeros(mesh.nboundelems)
    ref_vol = elementvolume_boundary(mesh.d)
    
    for el in eachindex(v)
        detJ = jacobian_boundary(mesh, el)
        v[el] = detJ * ref_vol
    end

    return v
end

"""
$(TYPEDSIGNATURES)

Returns barycenter vector of the d-dimensional reference element. 
"""
function elementbarycenter(d::Int64)
    v = zeros(Float64, d)
    v .= 1 / (d + 1)

    return v
end

"""
$(TYPEDSIGNATURES)

Returns barycenter of the given element in the given mesh. 
"""
function elementbarycenter(mesh::Mesh, element::Int64)
    J = base_jacobian(mesh, element)
    shift = mesh.Nodes[mesh.Elements[element][1]]
    return J * elementbarycenter(mesh.d) + shift 
end

"""
$(TYPEDSIGNATURES)

Returns vector of all barycenters of all element in the given mesh. 
"""
function elementbarycenter(mesh::Mesh)
    v = Array{Array{Float64,1},1}(undef, mesh.nelems)
    ref_bc = elementbarycenter(mesh.d)

    for el in eachindex(v)
        J = base_jacobian(mesh, el)
        shift = mesh.Nodes[mesh.Elements[el][1]]
        v[el] = J * ref_bc + shift
    end

    return v
end

"""
$(TYPEDSIGNATURES)
    
Returns diameter (i.e. longest edge length) of an element
spanned by nodes at the given coordinates.

The coordinates do not necessarily have to span an element of full dimension.
Could also be used for boundary elements.
"""
function elementdiameter(coords::Vector{Vector{Float64}})
    max = 0
    for i = length(coords):-1:2
        for j = i-1:-1:1
            val = norm(coords[i]-coords[j])
            if val > max
                max = val
            end
        end
    end

    return max
end

"""
$(TYPEDSIGNATURES)
    
Same as previous `$(FUNCTIONNAME)(...)`, but takes a mesh and set of nodes as arguments.
Extracts coordinates from the mesh and passes them to the base function.

Commonly the coordinates correspond to one (boundary) element in a mesh,
but not necessarily have to.
"""
function elementdiameter(mesh::Mesh, nodes::Vector{Int64})
    return elementdiameter(mesh.Nodes[nodes])
end

"""
$(TYPEDSIGNATURES)
    
Same as previous `$(FUNCTIONNAME)(...)`, but takes a mesh and an element id
as arguments. Extracts coordinates of the support nodes from the mesh
and passes them to the base function.
"""
function elementdiameter(mesh::Mesh, element::Int64)
    return elementdiameter(mesh, mesh.Elements[element])
end

"""
$(TYPEDSIGNATURES)
    
Similar ot previous `$(FUNCTIONNAME)(...)`, but returns vector of elementdiameters for all 
elements in the given mesh.
"""
function elementdiameter(mesh::Mesh)
    v = zeros(mesh.nelems)

    for el in eachindex(v)
        v[el] = elementdiameter(mesh, el)
    end

    return v
end

"""
$(TYPEDSIGNATURES)
    
Returns diameter (i.e. longest edge length) of a boundary element
spanned by nodes at the given coordinates.

Commonly the coordinates correspond to one boundary element in a mesh,
but not necessarily have to.

Since the number of nodes to span an element is not prescribed,
the functionality is identical to `elementdiameter(args)`.
However this function is kept for consistency and improved readability.
"""
function elementdiameter_boundary(coords::Vector{Vector{Float64}})
    return elementdiameter(coords)
end

"""
$(TYPEDSIGNATURES)
    
Same as previous `$(FUNCTIONNAME)(...)`, but takes a mesh and set of nodes as arguments.
Extracts coordinates from the mesh and passes them to the base function.

Commonly the coordinates correspond to one (boundary) element in a mesh,
but not necessarily have to.

Since the number of nodes to span an element is not prescribed,
the functionality is identical to `elementdiameter(args)`.
However this function is kept for consistency and improved readability.
"""
function elementdiameter_boundary(mesh::Mesh, nodes::Vector{Int64})
    return elementdiameter_boundary(mesh.Nodes[nodes])
end

"""
$(TYPEDSIGNATURES)
    
Same as previous `$(FUNCTIONNAME)(...)`, but takes a mesh and a boundary element id
as arguments. Extracts coordinates of the support nodes from the mesh
and passes them to the base function.
"""
function elementdiameter_boundary(mesh::Mesh, element::Int64)
    return elementdiameter_boundary(mesh, mesh.BoundaryElements[element])
end

"""
$(TYPEDSIGNATURES)
    
Similar ot previous `$(FUNCTIONNAME)(...)`, but returns vector of elementdiameters for all 
elements in the given mesh.
"""
function elementdiameter_boundary(mesh::Mesh)
    v = zeros(mesh.nboundelems)

    for el in eachindex(v)
        v[el] = elementdiameter_boundary(mesh, el)
    end

    return v
end

"""
$(TYPEDSIGNATURES)

Returns the radius of the circumscribed ball of the one-dimenional element
spanned by two one-dimensional coords. In one dimension, 
this reduces to half the elements length. 
"""
function circumscribedball1d(coords::Array{Array{Float64,1},1})
    return abs(coords[2][1] - coords[1][1]) / 2
end

"""
$(TYPEDSIGNATURES)

Returns the radius of the circumscribed ball of the two-dimenional element
spanned by three two-dimensional coords.
"""
function circumscribedball2d(coords::Array{Array{Float64,1},1})
    l12 = norm(coords[2] - coords[1])
    l13 = norm(coords[3] - coords[1])
    l23 = norm(coords[3] - coords[2])

    lc = 0.5 * (l12+l13+l23)
    area = sqrt(lc*(lc-l12)*(lc-l13)*(lc-l23))

    return l12*l13*l23 / (4*area)
end

"""
$(TYPEDSIGNATURES)

Returns the radius of the circumscribed ball of the three-dimenional element
spanned by four three-dimensional coords.
"""
function circumscribedball3d(coords::Array{Array{Float64,1},1})
    e12 = coords[2]-coords[1]
    e13 = coords[3]-coords[1]
    e14 = coords[4]-coords[1]
    
    J = [e12 e13 e14]
    v = [norm(e12)^2, norm(e13)^2, norm(e14)^2]

    return 0.5 * norm(inv(J)' * v)
end

"""
$(TYPEDSIGNATURES)

Returns the radius of the circumscribed ball for an element spanned by the given coordinates.
"""
function circumscribedball(coords::Array{Array{Float64,1},1})
    dim = length(coords)-1
    if dim == 1
        return circumscribedball1d(coords)
    elseif dim == 2
        return circumscribedball2d(coords)
    elseif dim == 3
        return circumscribedball3d(coords)
    else
        throw(ArgumentError("Unsuitable set of coordinates to span element."))
    end
end

"""
$(TYPEDSIGNATURES)
    
Same as previous `$(FUNCTIONNAME)(...)`, but takes a mesh and set of nodes as arguments.
Extracts coordinates from the mesh and passes them to the base function.

Commonly the coordinates correspond to one element in a mesh, but not necessarily have to.
"""
function circumscribedball(mesh::Mesh, nodes::Array{Int64,1})
    return circumscribedball(mesh.Nodes[nodes])
end

"""
$(TYPEDSIGNATURES)
    
Same as previous `$(FUNCTIONNAME)(...)`, but takes a mesh and an element id as arguments.
Extracts coordinates of the support nodes from the mesh
and passes them to the base function.
"""
function circumscribedball(mesh::Mesh, element::Int64)
    return circumscribedball(mesh, mesh.Elements[element])
end

"""
$(TYPEDSIGNATURES)

Returns a vector of the radii of the circumscribed balls of each element in the given mesh.
"""
function circumscribedball(mesh::Mesh)
    v = zeros(mesh.nelems)
    ref_vol = elementvolume(mesh.d)

    for el in eachindex(v)
        v[el] = circumscribedball(mesh, el)
    end

    return v
end

"""
$(TYPEDSIGNATURES)

Returns the radius of the inscribed ball in the one-dimenional element
spanned by two one-dimensional coords. In one dimension, 
this reduces to half the elements length. 
"""
function inscribedball1d(coords::Array{Array{Float64,1},1})
    return abs(coords[2][1] - coords[1][1]) / 2
end

"""
$(TYPEDSIGNATURES)

Returns the radius of the inscribed ball in the two-dimenional element
spanned by three two-dimensional coords.
"""
function inscribedball2d(coords::Array{Array{Float64,1},1})
    l12 = norm(coords[2] - coords[1])
    l13 = norm(coords[3] - coords[1])
    l23 = norm(coords[3] - coords[2])

    lc = 0.5 * (l12+l13+l23)
    area = sqrt(lc*(lc-l12)*(lc-l13)*(lc-l23))

    return area / lc
end

"""
$(TYPEDSIGNATURES)

Returns the radius of the inscribed ball in the three-dimenional element
spanned by four three-dimensional coords.
"""
function inscribedball3d(coords::Array{Array{Float64,1},1})
    e12 = coords[2]-coords[1]
    e13 = coords[3]-coords[1]
    e14 = coords[4]-coords[1]
    e23 = coords[3]-coords[2]
    e24 = coords[4]-coords[2]

    s123 = 0.5 * norm(cross(e12, e13))
    s124 = 0.5 * norm(cross(e12, e14))
    s134 = 0.5 * norm(cross(e13, e14))
    s234 = 0.5 * norm(cross(e23, e24))
    surface = s123 + s124 + s134 + s234

    J = [e12 e13 e14]
    volume =  det(J) * elementvolume(3)

    return 3 * volume / surface
end

"""
$(TYPEDSIGNATURES)

Returns the radius of the inscribed ball for an element spanned by the given coordinates.
"""
function inscribedball(coords::Array{Array{Float64,1},1})
    dim = length(coords)-1
    if dim == 1
        return inscribedball1d(coords)
    elseif dim == 2
        return inscribedball2d(coords)
    elseif dim == 3
        return inscribedball3d(coords)
    else
        throw(ArgumentError("Unsuitable set of coordinates to span element."))
    end
end

"""
$(TYPEDSIGNATURES)
    
Same as previous `$(FUNCTIONNAME)(...)`, but takes a mesh and set of nodes as arguments.
Extracts coordinates from the mesh and passes them to the base function.

Commonly the coordinates correspond to one element in a mesh, but not necessarily have to.
"""
function inscribedball(mesh::Mesh, nodes::Array{Int64,1})
    return inscribedball(mesh.Nodes[nodes])
end

"""
$(TYPEDSIGNATURES)
    
Same as previous `$(FUNCTIONNAME)(...)`, but takes a mesh and an element id as arguments.
Extracts coordinates of the support nodes from the mesh
and passes them to the base function.
"""
function inscribedball(mesh::Mesh, element::Int64)
    return inscribedball(mesh, mesh.Elements[element])
end

"""
$(TYPEDSIGNATURES)

Returns a vector of the radii of the inscribed balls of each element in the given mesh.
"""
function inscribedball(mesh::Mesh)
    v = zeros(mesh.nelems)
    ref_vol = elementvolume(mesh.d)

    for el in eachindex(v)
        v[el] = inscribedball(mesh, el)
    end

    return v
end

"""
$(TYPEDSIGNATURES)

Returns gridsize parameter h of the mesh,
which is given by the maximum of all element diameters.
Since we consider only tetrahedral elements,
the diameters correspond to the longest edge length.
Hence, the gridsize is the overall longest edge in the mesh.
"""
function gridsize(mesh::Mesh)
    return maximum(elementdiameter(mesh))
end

"""
$(TYPEDSIGNATURES)

Returns quasi-uniformity constant in the sense

```math
    \\min_{K \\in T_h} \\frac{h_K}{h} \\geq C
```

where ``h_K`` is the diameter of element ``K``
and ``h`` the gridsize of the mesh ``T_h``, i.e., the maximum over all diameters..
"""
function quasiuniformity(mesh::Mesh)
    hi = elementdiameter(mesh)
    h = maximum(hi)
    
    return minimum(hi ./ h)
end

"""
$(TYPEDSIGNATURES)

Returns shape regularity constant in the sense

```math
    \\max_{K \\in T_h} \\frac{h_K}{\\rho_K} \\leq C
```

where ``h_K`` is the diameter of element ``K``
and ``\\rho_K`` denotes the radius of the largest inscribed cirlce in ``h_K``.
"""
function shaperegularity(mesh::Mesh)
    hi = elementdiameter(mesh)
    ri = inscribedball(mesh)
    
    return maximum(hi ./ ri)
end

"""
$(TYPEDSIGNATURES)

Returns ratio of circumscribed to inscribed circle or sphere
for an element spanned by nodes at the given coordinates.
Also sometimes refered as aspectratio of the element.

Optimal ratios are equal to the respective dimension of the element.
Hence, the result can be normed with a factor 1/dim.
"""
function elementratio(coords::Array{Array{Float64,1},1})
    dim = length(coords)-1
    if dim == 1
        return 1
    elseif dim == 2
        return circumscribedball2d(coords) / inscribedball2d(coords)
    elseif dim == 3
        return circumscribedball3d(coords) / inscribedball3d(coords)
    else
        throw(ArgumentError("Unsuitable set of coordinates to span element."))
    end
end

"""
$(TYPEDSIGNATURES)

Same as previous `$(FUNCTIONNAME)(...)`,
but returns vector of ratios for all elements in the given mesh.
"""
function elementratio(mesh::Mesh)
    ratios = zeros(mesh.nelems)
    for el=1:mesh.nelems
        ratios[el] = elementratio(mesh.Nodes[mesh.Elements[el]])
    end

    return ratios
end

"""
$(TYPEDSIGNATURES)

Returns smallest interior angle for an element spanned by nodes at the given coordinates.
"""
function elementangle(coords::Array{Array{Float64,1},1})
    n = length(coords)
    if n == 2
        throw(ErrorException("Element angle for 1D currently not supported."))
    end
    
    min = pi
    for i = 1:n
        c = coords[i]
        for j = 1:n
            for k = 1:n
                if i != j && i != k && j != k
                    a = c - coords[j]
                    b = c - coords[k]

                    angle = acos(dot(a,b)/(norm(a)*norm(b)))
                    if angle < min
                        min = angle
                    end
                end
            end
        end
    end

    return min
end

"""
$(TYPEDSIGNATURES)

Same as previous `$(FUNCTIONNAME)(...)`,
but returns vector of angles for all elements in the given mesh.
"""
function elementangle(mesh::Mesh)
    angles = zeros(mesh.nelems)
    for el=1:mesh.nelems
        angles[el] = elementangle(mesh.Nodes[mesh.Elements[el]])
    end

    return angles 
end

"""
$(TYPEDSIGNATURES)
    
Returns the d-dimenional volume of the domain definded by the mesh. 
"""
function volume(mesh::Mesh)
    v = 0
    ref_vol = elementvolume(mesh.d)

    for el = 1:mesh.nelems
        detJ = det(base_jacobian(mesh, el))
        v += detJ * ref_vol
    end

    return v
end

"""
$(TYPEDSIGNATURES)
    
Returns the d-dimenional domain of the domain definded by the mesh. 
"""
function barycenter(mesh::Mesh)
    ref_vol = elementvolume(mesh.d)
    ref_bc = elementbarycenter(mesh.d)

    bc = zeros(Float64,mesh.d)
    vol = 0
    for el = 1:mesh.nelems
        J = base_jacobian(mesh, el)
        
        shift = mesh.Nodes[mesh.Elements[el][1]]
        ebc = J * ref_bc .+ shift
        evol = abs(det(J) * ref_vol)

        bc += evol .* ebc
        vol += evol
    end

    return bc ./ vol
end

"""
$(TYPEDSIGNATURES)

Returns width L of a strip that the meshed domain fits into.
"""
function stripwidth(mesh::Mesh)
    points = boundingbox(mesh)
    return maximum(points[2] .- points[1])
end

"""
$(TYPEDSIGNATURES)

Returns two nodes which span the bounding box of the mesh.
"""
function boundingbox(mesh::Mesh)
    min = Inf .* ones(mesh.d)
    max = -Inf .* ones(mesh.d)

    for v in mesh.Nodes
        for k = 1:mesh.d
            if v[k] < min[k]
                min[k] = v[k]
            elseif v[k] > max[k]
                max[k] = v[k]
            end
        end
    end
    return Array[min,max]
end
