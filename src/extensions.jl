"""
    barycenter(mesh::Mesh) -> Vector{Float64}
    
Returns vector of coordinates of the barycenter of the domain definded by the mesh. 
"""
function barycenter(mesh::Mesh)
    M = assemble_massmatrix(mesh)
    vol = volume(mesh)
    bc = zeros(Float64, mesh.d)
    for j = 1:mesh.d
        bc[j] =  sum(M * evaluate_mesh_function(mesh, x->x[j])) / vol
    end 
    
    return bc 
end

"""
    pnorm(p::Float64,v::AbstractVector{Float64}, mesh::Mesh; qdim::Int64=1, order::Int64=1) -> Float64

Returns \$L^p\$-norm of a FEM coefficient vector v over the nodes or quadrature points in the elements on the given mesh.
"""
function pnorm(p::Float64, v::AbstractVector{Float64}, mesh::Mesh; qdim::Int64=1, order::Int64=1)
    if p < 1
        throw(DomainError(p, "p has to be conjugated Hölder exponent 1 ≤ p ≤ ∞."))
    end
    
    if length(v) == mesh.nnodes * qdim
        if isinf(p)
            t = abs.(v)
            le = mesh.nnodes
        else
            E = assemble_basismatrix(mesh, qdim=qdim, order=order)
            t = abs.(E * v)
            le = mesh.nelems * length(quadrature_weights(mesh.d, order))
        end
    elseif mod(length(v), mesh.nelems*qdim) == 0
        t = abs.(v)

        nPoints = div(length(v), mesh.nelems * qdim)
        order = quadrature_order(mesh.d, nPoints)
        le = mesh.nelems * nPoints
    else
        throw(ArgumentError("The vector v does not have a valid length."))
    end

    if qdim != 1
        t = sqrt.(restrict_multivector(t.^2, le, qdim))
    end
    
    
    if isinf(p)       
        return maximum(t)
    else
        w = assemble_weightmultivector(mesh, qdim=1, order=order)

        if p == 1
            return sum(w .* t)
        else
            return sum(w .* t.^p).^(1/p)
        end
    end
end

"""
    qnorm(p::Float64,v::AbstractVector{Float64}, mesh::Mesh; qdim::Int64=1, order::Int64=1) -> Float64

Returns \$L^q\$-norm of a function f or its FEM coefficient vector v on the given mesh
with q being the conjugated exponent to p.
"""
function qnorm(p::Float64, v::AbstractVector{Float64}, mesh::Mesh; qdim::Int64=1, order::Int64=1)
    return pnorm(conjugated_exponent(p), v, mesh, qdim=qdim, order=order)
end


"""
    conjugated_exponent(p::Float64) -> Float64

Returns conjugated exponent q to p in the Hölder sense 1/p + 1/q = 1.
"""
function conjugated_exponent(p::Float64)
    if isinf(p)
        return 1.0
    elseif p == 1
        return Inf
    else
        return 1.0/(1-(1/p))
    end
end

"""
    pnorm_boundary(p::Float64, v::AbstractVector{Float64}, mesh::Mesh; boundaryElements::Set{Int64}=Set{Int64}(), qdim::Int64=1, order::Int64=1) -> Float64

Returns \$L^p\$-norm of a FEM coefficient vector v over the nodes or quadrature points in the boundary elements on the given mesh.
"""
function pnorm_boundary(p::Float64, v::AbstractVector{Float64}, mesh::Mesh; boundaryElements::Set{Int64}=Set{Int64}(), qdim::Int64=1, order::Int64=1)
    if p < 1
        throw(DomainError(p, "p has to be conjugated Hölder exponent 1 ≤ p ≤ ∞."))
    end

    if isempty(boundaryElements)
        boundaryElements = Set{Int64}(1 : mesh.nboundelems)
    end
    
    if length(v) == mesh.nnodes * qdim
        if isinf(p)
            t = abs.(v)
            le = mesh.nnodes
        else
            E = assemble_basismatrix_boundary(mesh, boundaryElements=boundaryElements, qdim=qdim, order=order)
            t = abs.(E * v)
            le = mesh.nboundelems * length(quadrature_weights(mesh.d-1, order))
        end
    elseif mod(length(v), mesh.nboundelems*qdim) == 0
        t = abs.(v)

        nPoints = div(length(v), mesh.nboundelems * qdim)
        order = quadrature_order(mesh.d-1, nPoints)
        le = mesh.nboundelems * nPoints
    end
    
    if qdim != 1
        t = sqrt.(restrict_multivector(t.^2, le, qdim))
    end

    if isinf(p)
        return maximum(t)
    else
        w = assemble_weightmultivector_boundary(mesh, qdim=1, order=order)

        if p == 1
            return sum(w.* t)
        else
            return sum(w .* t.^p).^(1/p)
        end
    end
end

"""
    qnorm_boundary(p::Float64, v::AbstractVector{Float64}, mesh::Mesh; boundaryElements::Set{Int64}=Set{Int64}(), qdim::Int64=1, order::Int64=1) -> Float64

Returns \$L^q\$-norm of a function f or its FEM coefficient vector v on the given mesh
with q being the conjugated exponent to p.
"""
function qnorm_boundary(p::Float64, v::AbstractVector{Float64}, mesh::Mesh; boundaryElements::Set{Int64}=Set{Int64}(), qdim::Int64=1, order::Int64=1)
    return pnorm_boundary(conjugated_exponent(p), v, mesh, boundaryElements=boundaryElements, qdim=qdim, order=order)
end
