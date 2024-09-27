"""
    pnorm(
        p::Float64,
        v::AbstractVector{Float64},
        mesh::Mesh; 
        qdim::Int64 = 1,
        order::Int64 = 1
    ) -> Float64

Returns \$L^p\$-norm of a FEM coefficient vector v over the nodes 
or quadrature points in the elements on the given mesh.
Note that for p=Inf, the values do not get inter- or extrapolated
and thus it might not be true maximum over the domain, but only over the given data.
"""
function pnorm(
    p::Float64,
    v::AbstractVector{Float64},
    mesh::Mesh; 
    qdim::Int64 = 1,
    order::Int64 = 1
)
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
    qnorm(
        p::Float64,
        v::AbstractVector{Float64},
        mesh::Mesh;
        qdim::Int64 = 1,
        order::Int64 = 1
    ) -> Float64

Returns \$L^q\$-norm of a FEM coefficient vector v over the nodes or quadrature points in
the elements on the given mesh with q being the conjugated exponent to p.
"""
function qnorm(
    p::Float64,
    v::AbstractVector{Float64},
    mesh::Mesh;
    qdim::Int64 = 1,
    order::Int64 = 1
)
    return pnorm(conjugated_exponent(p), v, mesh, qdim=qdim, order=order)
end


"""
$(TYPEDSIGNATURES)

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
    pnorm_boundary(
        p::Float64,
        v::AbstractVector{Float64},
        mesh::Mesh;
        boundaryElements::Set{Int64} = Set{Int64}(),
        qdim::Int64 = 1,
        order::Int64 = 1
    ) -> Float64

Returns \$L^p\$-norm of a FEM coefficient vector v over the nodes or quadrature points 
in the boundary elements on the given mesh.
Note that for p=Inf, the values do not get inter- or extrapolated
and thus it might not be true maximum over the domain, but only over the given data.
"""
function pnorm_boundary(
    p::Float64,
    v::AbstractVector{Float64},
    mesh::Mesh;
    boundaryElements::Set{Int64} = Set{Int64}(),
    qdim::Int64 = 1,
    order::Int64 = 1
)
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
            E = assemble_basismatrix_boundary(
                mesh,
                boundaryElements=boundaryElements,
                qdim=qdim,
                order=order
            )
            t = abs.(E * v)
            le = mesh.nboundelems * length(quadrature_weights_boundary(mesh.d, order))
        end
    elseif mod(length(v), mesh.nboundelems*qdim) == 0
        t = abs.(v)

        nPoints = div(length(v), mesh.nboundelems * qdim)
        order = quadrature_order(mesh.d-1, nPoints)
        le = mesh.nboundelems * nPoints
    else
        throw(ArgumentError("The vector v does not have a valid length."))
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
    qnorm_boundary(
        p::Float64,
        v::AbstractVector{Float64},
        mesh::Mesh;
        boundaryElements::Set{Int64} = Set{Int64}(),
        qdim::Int64 = 1,
        order::Int64 = 1
    ) -> Float64

Returns \$L^q\$-norm of a function f or its FEM coefficient vector v over the nodes or
quadrature points in the boundary elements on the given mesh
with q being the conjugated exponent to p.
"""
function qnorm_boundary(
    p::Float64,
    v::AbstractVector{Float64},
    mesh::Mesh;
    boundaryElements::Set{Int64} = Set{Int64}(),
    qdim::Int64 = 1,
    order::Int64 = 1
)
    return pnorm_boundary(
        conjugated_exponent(p),
        v,
        mesh,
        boundaryElements=boundaryElements,
        qdim=qdim,
        order=order
    )
end

"""
$(TYPEDSIGNATURES)

Implementation of the discrete \$L^2\$-norm based on a mass matrix.
Faster version of pnorm for p=2 by using that direct assembly of mass matrix
instead of basis matrices is possible.
Can be used for domain or boundary depending on the respective mass matrix.
"""
function twonorm(v::AbstractVector{Float64}, M::AbstractMatrix{Float64})
    return sqrt(v' * M * v)
end

"""
    twonorm(
        v::AbstractVector{Float64},
        mesh::Mesh; 
        qdim::Int64 = 1,
        order::Int64 = 3
    )

Same as previous $(FUNCTIONNAME)(...), but takes mesh and number of components as arguments.
Then assembles the mass matrix and passes it to the base function.
"""
function twonorm(
    v::AbstractVector{Float64},
    mesh::Mesh; 
    qdim::Int64 = 1,
    order::Int64 = 3
)
    M = assemble_massmatrix(mesh, qdim=qdim, order=order)

    return twonorm(v, M)
end

"""
    twonorm_boundary(
        v::AbstractVector{Float64},
        mesh::Mesh;
        boundaryElements::Set{Int64} = Set{Int64}(),
        qdim::Int64 = 1,
        order::Int64 = 3
    )

Implementation of the discrete \$L^2\$-norm on the boundary based on a mass matrix.
Faster version of pnorm for p=2 by using that direct assembly of mass matrix
instead of basis matrices is possible.
Assembles boundary mass matrix form mesh, boundary elements and the number of components.
Then passes it to the base implementation of `twonorm`, which can be used for boundary
terms as well.
"""
function twonorm_boundary(
    v::AbstractVector{Float64},
    mesh::Mesh;
    boundaryElements::Set{Int64} = Set{Int64}(),
    qdim::Int64 = 1,
    order::Int64 = 3
)
    M = assemble_massmatrix_boundary(
        mesh,
        boundaryElements = boundaryElements,
        qdim = qdim,
        order = order
    )

    return twonorm(v, M)
end
