"""
    gausslegendre_points(order::Int64) -> Array{Float64,1}
    
Returns coordinates of the Gauss-Legendre quadrature points on the default interval [-1,1] 
for exact integration of polynomials up to the order given by a non-negative integer.

The maximum order of accuary is 9. The function is only used internally and thus does not
perform an check on the given order, but will either provide 1st order for negative
arguments or 9th order if the argument is larger.
"""
function gausslegendre_points(order::Int64)
    if order <= 1
        return [0]
    elseif order <= 3
        val = 1/sqrt(3)
        return [-val, val]
    elseif order <= 5
        val = sqrt(3/5)
        return  [-val, 0, val]
    elseif order <= 7
        val = 2/7 * sqrt(6/5)
        val1 = sqrt(3/7 - val)
        val2 = sqrt(3/7 + val)
        return  [-val2, -val1, val1, val2]
    else
        val = 2 * sqrt(10/7)
        val1 = 1/3 * sqrt(5 - val)
        val2 = 1/3 * sqrt(5 + val)
        return  [-val2, -val1, 0, val1, val2]
    end
end

"""
    gausslegendre_weights(order::Int64) -> Array{Float64,1}
    
Returns weights of the Gauss-Legendre quadrature points on the default interval [-1,1] 
for exact integration of polynomials up to the ordergiven by a non-negative integer.

The maximum order of accuary is 9. The function is only used internally and thus does not
perform an check on the given order, but will either provide 1st order for negative
arguments or 9th order if the argument is larger.
"""
function gausslegendre_weights(order::Int64)
    if order <= 1
        return [2]
    elseif order <= 3
        return [1, 1]
    elseif order <= 5
        val = 5/9
        return  [val, 8/9, val]
    elseif order <= 7
        val = sqrt(30)
        val1 = (18 + val) / 36
        val2 = (18 - val) / 36
        return  [val2, val1, val1, val2]
    else
        val = 13 * sqrt(70)
        val1 = (322 + val) / 900
        val2 = (322 - val) / 900
        return  [val2, val1, 128/225, val1, val2]
    end
end

"""
    compute_coordinates_line(order::Int64) -> Array{Array{Float64,1},1}
    
Returns coordinates of the Gauss-Legendre quadrature points 
on the 1-dimensional FEM reference element [0,1] 
for exact integration of polynomials up to the given order.
"""
function compute_coordinates_line(order::Int64)
    p1 = gausslegendre_points(order)
    n = length(p1)

    c = Array{Array{Float64,1},1}(undef, n)

    r = 1
    for i = 1:n
        c[r] = [(1 + p1[i]) / 2]
        r = r + 1    
    end
    return c
end

"""
    compute_weights_line(order::Int64) -> Array{Float64,1}
    
Returns weights of the Gauss-Legendre quadrature points 
on the 1-dimensional FEM reference element [0,1] 
for exact integration of polynomials up to the given order.
"""
function compute_weights_line(order::Int64)    
    p1 = gausslegendre_points(order)
    w1 = gausslegendre_weights(order)
    n = length(p1)
    
    w = Array{Float64,1}(undef, n)

    r = 1
    for i = 1:n
        w[r] = 0.5 * w1[i]
        r = r + 1    
    end
    return w
end

"""
    compute_coordinates_triangle(order::Int64) -> Array{Array{Float64,1},1}
    
Returns coordinates of the Gauss-Legendre quadrature points 
on the 2-dimensional FEM reference element 
for exact integration of polynomials up to the given order.
"""
function compute_coordinates_triangle(order::Int64)
    order1 = order + 1
    p1 = gausslegendre_points(order1)
    n = length(p1)

    c = Array{Array{Float64,1},1}(undef, n*n)

    r = 1
    for i = 1:n
        for j = 1:n
            c[r] = [(1 + p1[i]) / 2 , (1 - p1[i]) * (1 + p1[j]) / 4]
            r = r + 1    
        end
    end
    return c
end

"""
    compute_weights_triangle(order::Int64) -> Array{Float64,1}
    
Returns weights of the Gauss-Legendre quadrature points 
on the 2-dimensional FEM reference element 
for exact integration of polynomials up to the given order.
"""
function compute_weights_triangle(order::Int64)
    order1 = order + 1    
    p1 = gausslegendre_points(order1)
    w1 = gausslegendre_weights(order1)
    n = length(p1)
    
    w = Array{Float64,1}(undef, n*n)

    r = 1
    for i = 1:n
        for j = 1:n
            w[r] = (1 - p1[i]) / 8 * w1[i] * w1[j]
            r = r + 1    
        end
    end
    return w
end

"""
    compute_coordinates_tetrahedron(order::Int64) -> Array{Array{Float64,1},1}
    
Returns coordinates of the Gauss-Legendre quadrature points 
on the 3-dimensional FEM reference element 
for exact integration of polynomials up to the given order.
"""
function compute_coordinates_tetrahedron(order::Int64)
    order1 = order + 2
    p1 = gausslegendre_points(order1)
    n = length(p1)
    
    c = Array{Array{Float64,1},1}(undef, n*n*n)
    
    r = 1
    for i = 1:n
        for j = 1:n
            for k = 1:n
                c[r] = [(1 + p1[i]) / 2,
                        (1 - p1[i]) * (1 + p1[j]) / 4,
                        (1 - p1[i]) * (1 - p1[j]) * (1 + p1[k]) / 8 ]
                r = r + 1
            end    
        end
    end
    return c
end

"""
    compute_weights_tetrahedron(order::Int64) -> Array{Float64,1}
    
Returns weights of the Gauss-Legendre quadrature points 
on the 3-dimensional FEM reference element 
for exact integration of polynomials up to the given order.
"""
function compute_weights_tetrahedron(order::Int64)
    order1 = order + 2
    p1 = gausslegendre_points(order1)
    w1 = gausslegendre_weights(order1)
    n = length(p1)

    w = Array{Float64,1}(undef, n*n*n)
    
    r = 1   
    for i = 1:n
        for j = 1:n
            for k = 1:n
                w[r] = (1 - p1[i])^2 * (1 - p1[j]) / 64 * w1[i] * w1[j] * w1[k]
                r = r + 1
            end    
        end
    end
    return w
end

"""
    parentcoordinates(x::Array{Float64,1}, id::Int64)

Maps coordinates given on the d-1 dimensional reference element to
the coordinates on the boundary of the d dimensional reference element
specified by the local boundary id.
This can be used to map boundary quadrature points to the respective parent coordinates.
"""
function parentcoordinates(x::Array{Float64,1}, id::Int64)
    dim = length(x) + 1
    
    if dim == 1
        if id == dim+1
            return [1.0]
        else
            return [0.0]
        end
    elseif dim == 2 
        if id == 3
            return [1-x[1], x[1]]
        elseif id == 2
            return [0.0, x[1]]
        else
            return [x[1], 0.0]
        end
    elseif dim == 3
        if id == 4
            return [1-x[1]-x[2], x[1], x[2]]
        elseif id == 3
            return [0.0, x[1], x[2]]
        elseif id == 2
            return [x[1], 0.0, x[2]]
        else
            return [x[1], x[2], 0.0]
        end
    else
        throw(DomainError(
            dim,
            "Unsupported dimension. Only dimensions 0,1,2 and 3 are supported."
        ))
    end
end

"""
    quadrature_points(d::Int64, order::Int64) -> Array{Array{Float64,1},1}
    
Returns coordinates of the Gauss-Legendre quadrature points 
on the d-dimensional FEM reference element 
for exact integration of polynomials up to the given order.
"""
function quadrature_points(d::Int64, order::Int64)
    if order < 0
        throw(ErrorException("Order $order not possible. Order has to be at least 0."))
    end
    if d < 0 || d > 3
        throw(ErrorException("Dimension $d not possible. " *
                                "Quadrature points are only available " *
                                "for 0D, 1D, 2D and 3D."))
    end

    if d == 0
        return [Array{Float64,1}()]    
    elseif d == 1
        if order <= 9
            return compute_coordinates_line(order)
        else
            throw(ErrorException("Order $order not possible. " *
                                    "Highest possible order for 1D is 9."))
        end
    elseif d == 2
        if order <= 1
            return [[1/3, 1/3]]
        elseif order <= 8
            return compute_coordinates_triangle(order)
        else
            throw(ErrorException("Order $order not possible. " *
                                    "Highest possible order for 2D is 8."))
        end
    elseif d == 3
        if order <= 1
            return [[1/4, 1/4, 1/4]]
        elseif order <= 7
            return compute_coordinates_tetrahedron(order)
        else
            throw(ErrorException("Order $order not possible. " *
                                    "Highest possible order for 3D is 7."))
        end
    end
end

"""
    quadrature_points(mesh::Mesh, order::Int64) -> Array{Array{Float64,1},1}
    
Returns global coordinates of the Gauss-Legendre quadrature points 
on each finite element in the given mesh
for exact integration of polynomials up to the given order.
"""
function quadrature_points(mesh::Mesh, order::Int64)
    quadX = quadrature_points(mesh.d, order)
    xle = length(quadX)
    Points = Array{Array{Float64,1},1}(undef,xle*mesh.nelems)
    
    for (el,nodes) in enumerate(mesh.Elements)
        shift = mesh.Nodes[nodes[1]]
        J = base_jacobian(mesh, nodes)
        
        for (i,r) in enumerate(quadX)
            Points[xle*(el-1)+i] = J * r + shift
        end
    end

    return Points
end

"""
    quadrature_points(element::Int64, mesh::Mesh, order::Int64) -> Array{Array{Float64,1},1}
    
Returns global coordinates of the Gauss-Legendre quadrature points 
on a specified element in the given mesh
for exact integration of polynomials up to the given order.
"""
function quadrature_points(element::Int64, mesh::Mesh, order::Int64)
    quadX = quadrature_points(mesh.d, order)
    Points = Array{Array{Float64,1},1}(undef, length(quadX))
    
    nodes = mesh.Elements[element]
    shift = mesh.Nodes[nodes[1]]
    J = base_jacobian(mesh, nodes)
        
    for (i,r) in enumerate(quadX)
        Points[i] = J * r + shift
    end

    return Points
end

"""
    quadrature_points_boundary(d::Int64, order::Int64) -> Array{Array{Float64,1},1}
    
Returns coordinates of the Gauss-Legendre quadrature points 
on the (d-1)-dimensional FEM reference element 
for exact integration of polynomials up to the given order.
"""
function quadrature_points_boundary(d::Int64, order::Int64)
    return quadrature_points(d-1,order)
end

"""
    quadrature_points_boundary(
        mesh::Mesh,
        order::Int64;
        boundaryElements::Set{Int64}=Set{Int64}()
    )-> Array{Array{Float64,1},1}
    
Returns global coordinates of the Gauss-Legendre quadrature points 
on each boundary element in the given mesh
for exact integration of polynomials up to the given order.
"""
function quadrature_points_boundary(
    mesh::Mesh,
    order::Int64;
    boundaryElements::Set{Int64}=Set{Int64}()
)
    elements = 1:mesh.nboundelems
    if !isempty(boundaryElements)
        elements = collect(boundaryElements)
    end

    quadX = quadrature_points_boundary(mesh.d, order)
    xle = length(quadX)
    Points = Array{Array{Float64,1},1}(undef,xle*mesh.nboundelems)

    for el in elements
        parentnodes = mesh.Elements[mesh.ParentElements[el]]
        parentboundary = mesh.ParentBoundaries[el]
        shift = mesh.Nodes[parentnodes[1]]
        J = base_jacobian(mesh, parentnodes)
        
        for (i,r) in enumerate(quadX)
            c = parentcoordinates(r, parentboundary)
            Points[xle*(el-1)+i] = J * c + shift
        end
    end
    return Points
end

"""
    quadrature_points_boundary(element::Int64, mesh::Mesh, order::Int64)
        -> Array{Array{Float64,1},1}
    
Returns global coordinates of the Gauss-Legendre quadrature points 
on a specified boundary element in the given mesh
for exact integration of polynomials up to the given order.
"""
function quadrature_points_boundary(element::Int64, mesh::Mesh, order::Int64)
    quadX = quadrature_points_boundary(mesh.d, order)
    Points = Array{Array{Float64,1},1}(undef,length(quadX))
    
    parentnodes = mesh.Elements[mesh.ParentElements[element]]
    parentboundary = mesh.ParentBoundaries[element]
    J = base_jacobian(mesh, parentnodes)
    shift = mesh.Nodes[parentnodes[1]]
        
    for (i,r) in enumerate(quadX)
        c = parentcoordinates(r, parentboundary)
        Points[i] = J * c + shift
    end

    return Points
end

"""
    quadrature_weights(d::Int64, order::Int64) -> Array{Float64,1}
    
Returns weights of the Gauss-Legendre quadrature points 
on the d-dimensional FEM reference element 
for exact integration of polynomials up to the given order.
"""
function quadrature_weights(d::Int64, order::Int64)
    if order < 0
        throw(ErrorException("Order $order not possible. Order has to be at least 0."))
    end
    if d < 0 || d > 3
        throw(ErrorException("Dimension $d not possible. " *
                                "Quadrature weights are only available " *
                                "for 1D, 2D and 3D."))
    end

    if d == 0
        return [1.0]    
    elseif d == 1
        if order <= 9
            return compute_weights_line(order)
        else
            throw(ErrorException("Order $order not possible. " *
                                    "Highest possible order for 1D is 9."))
        end
    elseif d == 2
        if order <= 1
            return [0.5]
        elseif order <= 8
            return compute_weights_triangle(order)
        else
            throw(ErrorException("Order $order not possible. " *
                                    "Highest possible order for 2D is 8."))
        end
    elseif d == 3
        if order <= 1
            return [1/6]
        elseif order <= 7
            return compute_weights_tetrahedron(order)
        else
            throw(ErrorException("Order $order not possible. " * 
                                    "Highest possible order for 3D is 7."))
        end
    end
end

"""
    quadrature_weights_boundary(d::Int64, order::Int64) -> Array{Float64,1}
    
Returns weights of the Gauss-Legendre quadrature points 
on the (d-1)-dimensional FEM reference element 
for exact integration of polynomials up to the given order.
"""
function quadrature_weights_boundary(d::Int64, order::Int64)
    return quadrature_weights(d-1, order)
end

"""
    quadrature_order(d::Int64, n::Int64) -> Int64
    
Returns highest quadrature order archived by using n points in d dimensions.
"""
function quadrature_order(d::Int64, n::Int64)
    return 2 * trunc(Int64, n^(1/d)) - 1
end

"""
    integral_over_reference_element(f::Function, d::Int64;order::Int64=1) -> Float64
    
Returns integral of function f over the d-dimensional FEM reference element 
computed with Gauss-Legendre quadrature exact at least for polynomials up to given order. 
"""
function integral_over_reference_element(f::Function, d::Int64; order::Int64=1)
    x = quadrature_points(d, order)
    w = quadrature_weights(d, order)

    return sum(w .* f.(x))
end

"""
    evaluate_quadrature_function(
        mesh::Mesh,
        f::Function;
        region::Set{Int64} = Set{Int64}(),
        order::Int64 = 1, 
        qdim::Int64 = 1
    )

Returns vector of function f evaluated at each elements quadrature points 
of the given mesh for the specified order.
"""
function evaluate_quadrature_function(
    mesh::Mesh,
    f::Function;
    region::Set{Int64} = Set{Int64}(),
    order::Int64 = 1, 
    qdim::Int64 = 1
)
    quadX = quadrature_points(mesh.d, order)
    xle = length(quadX)

    v = zeros(Float64, qdim * mesh.nelems * xle)

    elements = 1:mesh.nelems
    if !isempty(region)
        elements = collect(region)
    end

    if qdim == 1
        for i in elements
            coordinates = quadrature_points(i, mesh, order)
            for j = 1:xle
                id = xle*(i-1)+j
                v[id] = f(coordinates[j])
            end
        end
    else
        for i in elements
            coordinates = quadrature_points(i, mesh, order)
            for j = 1:xle
                s = xle*qdim*(i-1) + qdim*(j-1) + 1
                e = xle*qdim*(i-1) + qdim*j
                v[s:e] = f(coordinates[j])
            end
        end
    end

    return v
end

function evaluate_quadrature_function(
    mesh::Mesh,
    f::Function,
    region::Set{Domain};
    order::Int64 = 1, 
    qdim::Int64 = 1
) 
    evaluate_quadrature_function(
        mesh,
        f,
        region = extract_elements(region),
        order = order,
        qdim = qdim
    )
end

"""
    evaluate_quadrature_boundary(
        mesh::Mesh,
        f::Function;
        region::Set{Int64} = Set{Int64}(),
        order::Int64 = 1, 
        qdim::Int64 = 1 
    )
    evaluate_quadrature_function_boundary(
        mesh::Mesh,
        f::Function,
        region::Set{Boundary};
        order::Int64 = 1, 
        qdim::Int64 = 1
    )

Returns vector of function f evaluated at each boundary elements quadrature points 
of the given mesh for the specified order.
"""
function evaluate_quadrature_function_boundary(
    mesh::Mesh,
    f::Function;
    region::Set{Int64} = Set{Int64}(),
    order::Int64 = 1, 
    qdim::Int64 = 1
)
    quadX = quadrature_points_boundary(mesh.d, order)
    xle = length(quadX)

    v = zeros(Float64, qdim * mesh.nboundelems * xle)

    elements = 1:mesh.nboundelems
    if !isempty(region)
        elements = collect(region)
    end

    if qdim == 1
        for i in elements
            coordinates = quadrature_points_boundary(i, mesh, order)
            for j = 1:xle
                id = xle*(i-1)+j
                v[id] = f(coordinates[j])
            end
        end
    else
        for i in elements
            coordinates = quadrature_points_boundary(i, mesh, order)
            for j = 1:xle
                s = xle*qdim*(i-1) + qdim*(j-1) + 1
                e = xle*qdim*(i-1) + qdim*j 
                v[s:e] = f(coordinates[j])
            end
        end
    end

    return v
end

function evaluate_quadrature_function_boundary(
    mesh::Mesh,
    f::Function,
    region::Set{Boundary};
    order::Int64 = 1, 
    qdim::Int64 = 1
) 
    evaluate_quadrature_function_boundary(
        mesh,
        f,
        region = extract_elements(region),
        order = order,
        qdim = qdim
    )
end
