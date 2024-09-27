"""
$(TYPEDSIGNATURES)
    
Returns ii-th local basis function evaluated at x. 
Dimension is determined by length(x). 
"""
function phi(ii::Int64, x::AbstractVector)
    if ii == 1
        return 1.0 - sum(x)
    else
        return x[ii-1]
    end
end

"""
$(TYPEDSIGNATURES)
    
Gradient of ii-th local basis functions in d dimensions.
Note that the gradient is constant on an element.
"""
function grad_phi(d::Int64, ii::Int64)
    grad = zeros(d,1)
    if ii == 1
        grad .= -1.0
    else
        grad[ii-1] = 1.0
    end
    return grad
end

"""
    restrict_multivector(
        x::AbstractVector{Float64},
        m::Int64,
        qdim::Int64;
        block::Int64 = 1
    ) -> Vector{Float64}
    
Restricts a multivector of qdim×block×m elements for qdim components 
to the regular vector of m blocks of size block.
"""
function restrict_multivector(
    x::AbstractVector{Float64},
    m::Int64,
    qdim::Int64; 
    block::Int64 = 1
)
    if qdim == 1
        return copy(x)
    elseif block == 1
        v = zeros(m)
        for i=1:m
            v[i] = sum(x[qdim*(i-1)+1:qdim*i])
        end
        return v
    else
        v = zeros(block*m)
        for i = 1:m
            for j = 1:qdim
                v[block*(i-1) + j] = 
                    sum(x[(qdim*block*(i-1)+(j-1)*block+1):(qdim*block*(i-1)+block*j)])
            end
        end
        return v
    end
end 

"""
    prolong_multivector(
        x::AbstractVector{Float64},
        m::Int64,
        qdim::Int64; 
        block::Int64 = 1
    ) -> Vector{Float64}
    
Prolongates a vector of m blocks of size block to a multivector 
for qdim components of length qdim×block×m.
"""
function prolong_multivector(
    x::AbstractVector{Float64},
    m::Int64,
    qdim::Int64;
    block::Int64 = 1
)
    if qdim == 1
        return copy(x)
    elseif block == 1
        v = zeros(qdim*m)
        for i = 1:m
            v[qdim*(i-1)+1:qdim*i] .= x[i]
        end
        return v
    else
        v = zeros(qdim*block*m)
        for i = 1:m
            for j = 1:qdim
                v[(qdim*block*(i-1)+(j-1)*block+1):(qdim*block*(i-1)+block*j)] = 
                    x[(block*(i-1)+1):(block*i)]
            end
        end
        return v
    end
end

"""
    assemble_weightmultivector(
        mesh::Mesh;
        qdim::Int64 = 1,
        order::Int64 = 1
    ) -> Vector{Float64}
    
Returns the vector of weights for all elements of mesh and qdim components.
Weight is equal to volume if 1st order integration by mid-point rule is used.
"""
function assemble_weightmultivector(
    mesh::Mesh;
    qdim::Int64 = 1,
    order::Int64 = 1
)
    quadW = quadrature_weights(mesh.d, order)
    wle = length(quadW)

    m = mesh.nelems
    omega = zeros(Float64, m * wle)
    for i = 1:m
        (detJ, _) = jacobian(mesh, mesh.Elements[i])
        omega[wle*(i-1)+1:wle*i] += abs(detJ) .* quadW
    end
    
    return prolong_multivector(omega, m*wle, qdim)
end

"""
    assemble_weightmultivector_boundary(
        mesh::Mesh;
        qdim::Int64 = 1,
        order::Int64 = 1
    ) -> Vector{Float64}
    
Returns the vector of weights for all boundary elements of mesh and qdim components.
Weight is equal to volume if 1st order integration by mid-point rule is used.
"""
function assemble_weightmultivector_boundary(
    mesh::Mesh;
    qdim::Int64 = 1,
    order::Int64 = 1
)
    quadW = quadrature_weights_boundary(mesh.d, order)
    wle = length(quadW)

    m = mesh.nboundelems
    omega = zeros(Float64, m*wle)
    for i = 1:m
        detJ = jacobian_boundary(mesh, mesh.BoundaryElements[i])
        omega[wle*(i-1)+1:wle*i] += abs(detJ) .* quadW
    end

    return prolong_multivector(omega, m*wle, qdim)
end

"""
    assemble_derivativematrix(mesh::Mesh; qdim::Int64=1) -> SparseMatrixCSC{Float64, Int64}
    
Returns the discrete derivative matrix for all elements of mesh and qdim components.
"""
function assemble_derivativematrix(mesh::Mesh; qdim::Int64=1)   
    AA = zeros(Float64, mesh.d * qdim * mesh.nelems * (mesh.d+1))
    II = zeros(Int64, length(AA))
    JJ = zeros(Int64, length(AA))
    n = 0

    for el = 1:mesh.nelems
        nodes = mesh.Elements[el]
        (_, J) = jacobian(mesh, nodes)

        elemMat = zeros(Float64, mesh.d, mesh.d+1)
        for i = 1:(mesh.d+1)
            gradi = J * grad_phi(mesh.d,i)
            for k = 1:mesh.d
                elemMat[k,i] += gradi[k]
            end
        end
  
        for k = 1:mesh.d
            for i = 1:(mesh.d+1)
                for r = 1:qdim
                    n = n + 1
                    II[n] = mesh.d * qdim * (el-1) + qdim * (k-1) + r
                    JJ[n] = qdim * (nodes[i]-1) + r
                    AA[n] = elemMat[k,i]
                end
            end
        end
    end
    D = sparse(II[1:n], JJ[1:n], AA[1:n])
    return dropzeros!(D)
end

"""
    assemble_laplacian(
        mesh::Mesh;
        qdim::Int64 = 1
    ) -> SparseMatrixCSC{Float64, Int64}

Returns the Laplacian stiffness matrix for all elements of mesh and qdim components.
"""
function assemble_laplacian(
    mesh::Mesh;
    qdim::Int64 = 1
)
    AA = zeros(Float64, qdim^2 * mesh.nelems * (mesh.d+1)^2)
    II = zeros(Int64, length(AA))
    JJ = zeros(Int64, length(AA))
    n = 0

    for el = 1:mesh.nelems
        nodes = mesh.Elements[el]
        (detJ, J) = jacobian(mesh, nodes)
        
        elemMat = zeros(mesh.d+1,mesh.d+1)
        for i = 1:(mesh.d+1)
            gradi = J * grad_phi(mesh.d,i)
            for j = 1:(mesh.d+1)
                gradj = J * grad_phi(mesh.d,j)
                for k = 1:mesh.d
                    elemMat[i,j] += gradi[k] * gradj[k]
                end
            end
        end
        elemMat *= detJ / factorial(mesh.d) 

        for i = 1:(mesh.d+1)
            for j = 1:(mesh.d+1)
                for r = 1:qdim
                    n = n+1
                    II[n] = qdim * (nodes[i]-1) + r
                    JJ[n] = qdim * (nodes[j]-1) + r
                    AA[n] = elemMat[i,j]
                end
            end
        end
    end
    L = sparse(II[1:n],JJ[1:n],AA[1:n])
    return dropzeros!(L)
end

"""
    assemble_laplacian(
        D::SparseMatrixCSC{Float64, Int64},
        w::AbstractVector{Float64}
    ) -> SparseMatrixCSC{Float64, Int64}

Returns the Laplacian stiffness matrix for all elements of a mesh and qdim components.

Takes existing derivative matrix D and weight vector w as arguments.
Can yield performance benefits compared to direct assembly if D already exists.
Note that number of components and quadrature order in D and w have to match.
"""
function assemble_laplacian(
    D::SparseMatrixCSC{Float64, Int64},
    w::AbstractVector{Float64}
)
    L = D' * Diagonal(w) * D
    dropzeros!(L)
    return L
end

"""
    assemble_derivativematrix_boundary(
        mesh::Mesh; 
        boundaryElements::Set{Int64} = Set{Int64}(), 
        qdim::Int64 = 1
    ) -> SparseMatrixCSC{Float64, Int64}
    
Returns the discrete derivative matrix for all or specied boundary elements 
of mesh and qdim components.

The implementation is based on the derivate tensor for the corresponding full element
and the fact that the gradient ist constant on the element for linear elements. 
"""
function assemble_derivativematrix_boundary(
    mesh::Mesh; 
    boundaryElements::Set{Int64} = Set{Int64}(),
    qdim::Int64 = 1
)
    if isempty(boundaryElements)
        boundaryElements = Set{Int64}(1 : mesh.nboundelems)
    end

    AA = zeros(Float64, mesh.d * qdim * mesh.nboundelems * (mesh.d+1))
    II = zeros(Int64, length(AA))
    JJ = zeros(Int64, length(AA))
    n = 0

    for bel in boundaryElements
        nodes = mesh.Elements[mesh.ParentElements[bel]]
        (_, J) = jacobian(mesh, nodes)

        elemMat = zeros(Float64, mesh.d, mesh.d+1)
        for i = 1:(mesh.d+1)
            gradi = J * grad_phi(mesh.d,i)
            for k = 1:mesh.d
                elemMat[k,i] += gradi[k]
            end
        end

        for k = 1:mesh.d
            for i = 1:(mesh.d+1)
                for r = 1:qdim
                    n = n + 1
                    II[n] = mesh.d * qdim * (bel-1) + qdim * (k-1) + r
                    JJ[n] = qdim * (nodes[i]-1) + r
                    AA[n] = elemMat[k,i]
                end
            end
        end
    end
    
    D = sparse(II[1:n], JJ[1:n], AA[1:n], mesh.d*qdim*mesh.nboundelems, qdim*mesh.nnodes)
    return dropzeros!(D)
end

"""
    assemble_normalderivativematrix(
        mesh::Mesh;
        boundaryElements::Set{Int64} = Set{Int64}(),
        qdim::Int64 = 1
    ) -> SparseMatrixCSC{Float64, Int64}
    
Returns the discrete normal derivative matrix for all or specied boundary elements 
of mesh and qdim components.

The implementation is based on the derivate tensor for the corresponding full element
and the fact that the gradient ist constant on the element for linear elements. 
"""
function assemble_normalderivativematrix(
    mesh::Mesh; 
    boundaryElements::Set{Int64} = Set{Int64}(),
    qdim::Int64 = 1
)
    if isempty(boundaryElements)
        boundaryElements = Set{Int64}(1 : mesh.nboundelems)
    end
    AA = zeros(Float64, qdim * mesh.nboundelems * (mesh.d+1))
    II = zeros(Int64, length(AA))
    JJ = zeros(Int64, length(AA))
    n = 0

    for bel in boundaryElements
        nodes = mesh.Elements[mesh.ParentElements[bel]]
        (_, J) = jacobian(mesh, nodes)
        eta = outernormalvector(mesh, bel, J)

        elemVec = zeros(Float64, mesh.d+1)
        for i = 1:(mesh.d+1)
            gradi = J * grad_phi(mesh.d,i)
            for k = 1:mesh.d
                elemVec[i] += gradi[k] * eta[k]
            end
        end

        for i = 1:(mesh.d+1)
            for r = 1:qdim
                n = n + 1
                II[n] = qdim * (bel-1) + r
                JJ[n] = qdim * (nodes[i]-1) + r
                AA[n] = elemVec[i]
            end
        end
    end
    
    D = sparse(II[1:n], JJ[1:n], AA[1:n], qdim*mesh.nboundelems, qdim*mesh.nnodes)
    return dropzeros!(D)
end

"""
    assemble_basismatrix(
        mesh::Mesh;
        qdim::Int64 = 1,
        order::Int64 = 3
    ) -> SparseMatrixCSC{Float64, Int64}
    
Returns the discrete basis matrix with given local integration order 
for all elements of mesh and qdim components.
"""
function assemble_basismatrix(
    mesh::Mesh;
    qdim::Int64 = 1,
    order::Int64 = 3
)       
    quadX = quadrature_points(mesh.d, order)
    xle = length(quadX)
    
    AA = zeros(Float64, xle * qdim * mesh.nelems * (mesh.d+1))
    II = zeros(Int64, length(AA))
    JJ = zeros(Int64, length(AA))
    n = 0
    
    for el = 1:mesh.nelems
        nodes = mesh.Elements[el]

        elemMat = zeros(Float64, xle, mesh.d+1)
        for l = 1:xle    
            for i = 1:(mesh.d+1)
                elemMat[l,i] += phi(i, quadX[l])
            end
        end

        for l = 1:xle
            for i = 1:(mesh.d+1)
                for r = 1:qdim
                    n = n + 1
                    II[n] = xle * qdim * (el-1) + qdim * (l-1) + r
                    JJ[n] = qdim * (nodes[i]-1) + r
                    AA[n] = elemMat[l,i]
                end
            end
        end
    end
    E = sparse(II[1:n], JJ[1:n], AA[1:n])
    dropzeros!(E)
    return E
end

"""
    assemble_massmatrix(
        mesh::Mesh;
        qdim::Int64 = 1,
        order::Int64 = 1
    ) -> SparseMatrixCSC{Float64, Int64}

Returns the mass matrix with given local integration order for all elements 
of mesh and qdim components.
"""
function assemble_massmatrix(
    mesh::Mesh;
    qdim::Int64 = 1,
    order::Int64 = 3
)
    AA = zeros(Float64, qdim^2 * mesh.nelems * (mesh.d+1)^2)
    II = zeros(Int64, length(AA))
    JJ = zeros(Int64, length(AA))
    n = 0

    quadX = quadrature_points(mesh.d, order)
    quadW = quadrature_weights(mesh.d, order)

    for el = 1:mesh.nelems
        nodes = mesh.Elements[el]
        (detJ, _) = jacobian(mesh, nodes)
        elemMat = zeros(mesh.d+1, mesh.d+1)
        for i = 1:(mesh.d+1)
            for j = 1:(mesh.d+1)
                for (q, x) in enumerate(quadX)
                    elemMat[i,j] += phi(i, x) * phi(j, x) * quadW[q] * detJ
                end
            end
        end

        for d = 1:qdim
            for i = 1:(mesh.d+1)
                for j = 1:(mesh.d+1)
                    n = n+1
                    II[n] = qdim*(nodes[i]-1)+d
                    JJ[n] = qdim*(nodes[j]-1)+d
                    AA[n] = elemMat[i,j]
                end
            end
        end
    end

    M = sparse(II[1:n],JJ[1:n],AA[1:n])
    dropzeros!(M)
    return M
end

"""
$(TYPEDSIGNATURES)

Returns the mass matrix with given local integration order for all elements 
of mesh and qdim components.

Takes existing basis matrix E and weight vector w as arguments.
Can yield performance benefits compared to direct assembly if E already exists.
Note that number of components and quadrature order in E and w have to match.
"""
function assemble_massmatrix(
    E::SparseMatrixCSC{Float64, Int64}, 
    w::AbstractVector{Float64}
)
    M = E' * Diagonal(w) * E
    dropzeros!(M)
    return M
end


"""
    assemble_basismatrix_boundary(
        mesh::Mesh; 
        boundaryElements::Set{Int64} = Set{Int64}(),
        qdim::Int64 = 1,
        order::Int64 = 1
    ) -> SparseMatrixCSC{Float64, Int64}
    
Returns the discrete basis matrix with given local integration order
for all or specified boundary elements of mesh and qdim components.
"""
function assemble_basismatrix_boundary(
    mesh::Mesh; 
    boundaryElements::Set{Int64} = Set{Int64}(), 
    qdim::Int64 = 1,
    order::Int64 = 3
)
    if isempty(boundaryElements)
        boundaryElements = Set{Int64}(1 : mesh.nboundelems)
    end

    quadX = quadrature_points_boundary(mesh.d, order)
    xle = length(quadX)
    
    AA = zeros(Float64, xle * qdim * mesh.nboundelems * mesh.d)
    II = zeros(Int64, length(AA))
    JJ = zeros(Int64, length(AA))
    n = 0

    for bel in boundaryElements
        nodes = mesh.BoundaryElements[bel]

        elemMat = zeros(Float64, xle, mesh.d)
        for l = 1:xle    
            for i = 1:mesh.d
                elemMat[l,i] += phi(i, quadX[l])
            end
        end

        for l = 1:xle
            for i = 1:mesh.d
                for r = 1:qdim
                    n = n + 1
                    II[n] = xle * qdim * (bel-1) + qdim * (l-1) + r
                    JJ[n] = qdim * (nodes[i]-1) + r
                    AA[n] = elemMat[l,i]
                end
            end
        end
    end
    E = sparse(II[1:n], JJ[1:n], AA[1:n], xle*qdim*mesh.nboundelems, qdim*mesh.nnodes)
    dropzeros!(E)
    return E
end

"""
    assemble_massmatrix_boundary(
        mesh::Mesh; 
        boundaryElements::Set{Int64} = Set{Int64}(), 
        qdim::Int64 = 1,
        order::Int64 = 1
    ) -> SparseMatrixCSC{Float64, Int64}

Returns the mass matrix with given local integration order 
for all or specified boundary elements of mesh and image dimension qdim components.
"""
function assemble_massmatrix_boundary(
    mesh::Mesh; 
    boundaryElements::Set{Int64} = Set{Int64}(), 
    qdim::Int64 = 1,
    order::Int64 = 3
)
    if isempty(boundaryElements)
        boundaryElements = Set{Int64}(1 : mesh.nboundelems)
    end

    AA = zeros(Float64, qdim^2 * mesh.nelems * (mesh.d)^2)
    II = zeros(Int64, length(AA))
    JJ = zeros(Int64, length(AA))
    n = 0

    quadX = quadrature_points_boundary(mesh.d, order)
    quadW = quadrature_weights_boundary(mesh.d, order)

    for bel in boundaryElements
        nodes = mesh.BoundaryElements[bel]
        detJ = jacobian_boundary(mesh, nodes)
        elemMat = zeros(mesh.d, mesh.d)
        for i = 1:mesh.d
            for j = 1:mesh.d
                for (q, x) in enumerate(quadX)
                    elemMat[i,j] += phi(i, x) * phi(j, x) * quadW[q] * detJ
                end
            end
        end

        for d = 1:qdim
            for i = 1:mesh.d
                for j = 1:mesh.d
                    n = n+1
                    II[n] = qdim*(nodes[i]-1)+d
                    JJ[n] = qdim*(nodes[j]-1)+d
                    AA[n] = elemMat[i,j]
                end
            end
        end
    end

    return sparse(II[1:n],JJ[1:n],AA[1:n], qdim*mesh.nnodes, qdim*mesh.nnodes)
end

"""
$(TYPEDSIGNATURES)

Returns the mass matrix with given local integration order 
for all or specified boundary elements of mesh and image dimension qdim components.

Takes existing boundary basis matrix E and weight vector w as arguments.
Can yield performance benefits compared to direct assembly if E already exists.
Note that number of components and quadrature order in E and w have to match.
"""
function assemble_massmatrix_boundary(
    E::SparseMatrixCSC{Float64, Int64},
    w::AbstractVector{Float64}
)
    M = E' * Diagonal(w) * E
    dropzeros!(M)
    return M
end

"""
    assemble_cubicterm(
        mesh::Mesh,
        y::AbstractVector;
        order::Int64 = 3
    ) -> SparseMatrixCSC{Float64, Int64}

Returns the cubic term of the standard semilinear parabolic equation.
"""
function assemble_cubicterm(
    mesh::Mesh,
    y::AbstractVector;
    order::Int64 = 3
)
    V = zeros(mesh.nnodes)

    quadX = quadrature_points(mesh.d, order)
    quadW = quadrature_weights(mesh.d, order)
    
    for el = 1:mesh.nelems
        nodes = mesh.Elements[el]
        (detJ, _) = jacobian(mesh, nodes)

        y_cubic = zeros(length(quadW))
        for (q, x) in enumerate(quadX)
            y_quad = 0
            for i = 1:mesh.d+1
                y_quad += y[nodes[i]] * phi(i, x)
            end
            y_cubic[q] = y_quad^3
        end

        for i = 1:mesh.d+1
            for (q, x) in enumerate(quadX)
                V[nodes[i]] += y_cubic[q] * phi(i, x) * quadW[q] * detJ
            end
        end
    end

    return V
end

"""
    assemble_cubicderivativematrix(
        mesh::Mesh,
        y::AbstractVector;
        order::Int64 = 3
    ) -> SparseMatrixCSC{Float64, Int64}

Returns the linearization of the cubic term of the standard semilinear elliptic equation.
"""
function assemble_cubicderivativematrix(
    mesh::Mesh,
    y::AbstractVector;
    order::Int64 = 3
)
    AA = zeros(Float64, mesh.nelems * (mesh.d+1)^2)
    II = zeros(Int64, length(AA))
    JJ = zeros(Int64, length(AA))
    n = 0

    quadX = quadrature_points(mesh.d, order)
    quadW = quadrature_weights(mesh.d, order)

    for el = 1:mesh.nelems
        nodes = mesh.Elements[el]
        (detJ, J) = jacobian(mesh, nodes)

        y_quadratic = zeros(length(quadW))
        for (q, x) in enumerate(quadX)
            y_quad = 0
            for i = 1:mesh.d+1
                y_quad += y[nodes[i]] * phi(i, x)
            end
            y_quadratic[q] = y_quad^2
        end

        elemMat = zeros(mesh.d+1,mesh.d+1)
        for i = 1:mesh.d+1
            for j = 1:mesh.d+1
                for (q, x) in enumerate(quadX)
                    elemMat[i,j] += 3.0 * y_quadratic[q] * 
                                    phi(i, x) * phi(j, x) * quadW[q] * detJ
                end
            end
        end

        for i = 1:mesh.d+1
            for j = 1:mesh.d+1
                n = n+1
                II[n] = nodes[i]
                JJ[n] = nodes[j]
                AA[n] = elemMat[i,j]
            end
        end
    end

    return sparse(II[1:n],JJ[1:n],AA[1:n])
end

"""
    assemble_cubicsecondderivativematrix(
        mesh::Mesh,
        y::AbstractVector,
        p::AbstractVector;
        order::Int64 = 3
    ) -> SparseMatrixCSC{Float64, Int64}

Returns the second derivative of the cubic term 
of the standard semilinear elliptic equation around the state y.
"""
function assemble_cubicsecondderivativematrix(
    mesh::Mesh,
    y::AbstractVector, 
    p::AbstractVector;
    order::Int64 = 3
)
    AA = zeros(Float64, mesh.nelems * (mesh.d+1)^2)
    II = zeros(Int64, length(AA))
    JJ = zeros(Int64, length(AA))
    n = 0

    quadX = quadrature_points(mesh.d, order)
    quadW = quadrature_weights(mesh.d, order)

    for el = 1:mesh.nelems
        nodes = mesh.Elements[el]
        (detJ, J) = jacobian(mesh, nodes)
        
        y_quad = zeros(length(quadW))
        p_quad = zeros(length(quadW))
        for i = 1:3
            for (q, x) in enumerate(quadX)
                y_quad[q] += y[nodes[i]] * phi(i, x)
                p_quad[q] += p[nodes[i]] * phi(i, x)
            end
        end

        elemMat = zeros(mesh.d+1, mesh.d+1)
        for i = 1:mesh.d+1
            for j = 1:mesh.d+1
                for (q, x) in enumerate(quadX)
                    elemMat[i,j] += 6.0 * y_quad[q] * p_quad[q] *
                                    phi(i, x) * phi(j, x) * quadW[q] * detJ
                end
            end
        end

        for i = 1:mesh.d+1
            for j = 1:mesh.d+1
                n = n+1
                II[n] = nodes[i]
                JJ[n] = nodes[j]
                AA[n] = elemMat[i,j]
            end
        end
    end

    return sparse(II[1:n],JJ[1:n],AA[1:n])
end

"""
$(TYPEDSIGNATURES)

Returns the local strain tensor for linear elasticity with gradient as argument.
"""
function strain_tensor(
    grad::AbstractMatrix
) :: AbstractMatrix
    return 0.5 * (grad + grad')
end

"""
$(TYPEDSIGNATURES)

Returns the local stress tensor for linear elasticity 
with constant coefficients λ, μ and gradient as argument.
"""
function stress_tensor(
    grad::AbstractMatrix,
    lambda::Float64,
    mu::Float64
) :: AbstractMatrix
    d = size(grad, 1)
    epsilon = strain_tensor(grad)
    return lambda * tr(epsilon) * Matrix{Float64}(I, d, d) + 2.0 * mu * epsilon
end

"""
$(TYPEDSIGNATURES)

Returns the linear elasticity stiffness matrix with constant coefficients λ and μ 
for all elements of mesh and image dimension qdim.
"""
function assemble_elasticity(
    mesh::Mesh,
    lambda::Float64,
    mu::Float64
)
    qdim = mesh.d

    AA = zeros(Float64, mesh.nelems * qdim^2 * (mesh.d+1)^2)
    II = zeros(Int64, length(AA))
    JJ = zeros(Int64, length(AA))
    n = 0

    for el = 1:mesh.nelems
        nodes = mesh.Elements[el]
        (detJ, J) = jacobian(mesh, nodes)
        elemMat = zeros(qdim*(mesh.d+1), qdim*(mesh.d+1))
        for i = 1:(mesh.d+1)
            for j = 1:(mesh.d+1)
                for ic = 1:qdim
                    for jc = 1:qdim
                        grad_i = zeros(mesh.d, mesh.d)
                        grad_i[:,ic] = J * grad_phi(mesh.d, i)
                        grad_j = zeros(mesh.d, mesh.d)
                        grad_j[:,jc] = J * grad_phi(mesh.d, j)
                        elemMat[qdim*(i-1)+ic, qdim*(j-1)+jc] += 
                            dot(stress_tensor(grad_j, lambda, mu), strain_tensor(grad_i))
                    end
                end
            end
        end
        elemMat *= detJ / factorial(mesh.d)

        for i = 1:(mesh.d+1)
            for j = 1:(mesh.d+1)
                for ic = 1:qdim
                    for jc = 1:qdim
                        n = n+1
                        II[n] = qdim * (nodes[i]-1) + ic
                        JJ[n] = qdim * (nodes[j]-1) + jc
                        AA[n] = elemMat[qdim*(i-1)+ic, qdim*(j-1)+jc]
                    end
                end
            end
        end
    end

    return sparse(II[1:n],JJ[1:n],AA[1:n])
end

"""
    assemble_dirichletcondition!(
        A::SparseMatrixCSC{Float64, Int64},
        DI::Set{Int64}; 
        rhs = [],
        bc = [],
        qdim::Int64 = 1,
        insert = 1.0
    )

Modify a stiffness matrix and a right hand side according to the given Dirichlet conditions.

DI has to be the set of node indices for which the condition should be active.
For vector valued states either DI can be set to each component that should have a
Dirichlet condtion or qdim is set, if all components should have the condition.
If the rhs shall be updated bc needs to be specified explicitly as well.
The value insert is put as diagonal element. Usually you want a 1.0 here.
"""
function assemble_dirichletcondition!(
    A::SparseMatrixCSC{Float64, Int64},
    DI::Set{Int64}; 
    rhs = [],
    bc = [],
    qdim::Int64 = 1,
    insert = 1.0
)
    if rhs != [] && bc != []
        for i in DI
            for d = 1:qdim
                ii = qdim * (i-1) + d
                bcind = A[:, ii].nzind
                rhs[bcind] -= A[bcind, ii] * bc[ii]
            end
        end
        for i in DI
            for d = 1:qdim
                ii = qdim * (i-1) + d
                rhs[ii] = bc[ii]
            end
        end
    end

    for i in DI
        for d = 1:qdim
            ii = qdim * (i-1) + d
            A[ii, A[ii,:].nzind] *= 0.0;
            A[:, ii] *= 0.0;
            A[ii, ii] = insert;
        end
    end
    dropzeros!(A)
end

"""
    assemble_dirichletcondition!(
        A::SparseMatrixCSC{Float64, Int64},
        DI::Set{Boundary};
        rhs = [],
        bc = [],
        qdim::Int64 = 1,
        insert = 1.0
    )

Same as previous `$(FUNCTIONNAME)(...)`, but takes `Set{Boundary}` instead of `Set{Int64}`
as argument for the Dirichlet nodes.
Thus extracts the nodes first and then passes them to the base function.
"""
function assemble_dirichletcondition!(
    A::SparseMatrixCSC{Float64, Int64},
    DI::Set{Boundary};
    rhs = [],
    bc = [],
    qdim::Int64 = 1,
    insert = 1.0
)
    assemble_dirichletcondition!(
        A,
        extract_nodes(DI);
        rhs = rhs,
        bc = bc,
        qdim = qdim,
        insert = insert
    )
end

"""
    assemble_dirichletcondition_rhs!(
        A::SparseMatrixCSC{Float64, Int64},
        rhs::AbstractVector{Float64},
        DI::Set{Int64},
        bc::AbstractVector{Float64};
        qdim::Int64 = 1
    )

Modify a right hand side according to the given Dirichlet conditions.
Behaviour is similar to `assemble_dirichletcondition!(...)` however the system matrix A is
not updated. Can be relevant for iterative algorithm, where the system matrix is constant
and only the right hand side changes. Then one can store the modified matrix and
only assemble the right hand side in every iteration.

DI has to be the set of node indices for which the condition should be active.
For vector valued states either DI can be set to each component that should have a
Dirichlet condtion or qdim is set, if all components should have the condition.
"""
function assemble_dirichletcondition_rhs!(
    A::SparseMatrixCSC{Float64, Int64},
    rhs::AbstractVector{Float64},
    DI::Set{Int64},
    bc::AbstractVector{Float64};
    qdim::Int64 = 1
)
    if rhs != [] && bc != []
        for i in DI
            for d = 1:qdim
                ii = qdim * (i-1) + d
                bcind = A[:, ii].nzind
                rhs[bcind] -= A[bcind, ii] * bc[ii]
            end
        end
        for i in DI
            for d = 1:qdim
                ii = qdim * (i-1) + d
                rhs[ii] = bc[ii]
            end
        end
    end
end

"""
    assemble_dirichletcondition_rhs!(
        A::SparseMatrixCSC{Float64, Int64},
        rhs::AbstractVector{Float64},
        DI::Set{Boundary},
        bc::AbstractVector{Float64};
        qdim::Int64 = 1
    )

Same as previous `$(FUNCTIONNAME)(...)`, but takes `Set{Boundary}` instead of `Set{Int64}`
as argument for the Dirichlet nodes.
Thus extracts the nodes first and then passes them to the base function.
"""
function assemble_dirichletcondition_rhs!(
    A::SparseMatrixCSC{Float64, Int64},
    rhs::AbstractVector{Float64},
    DI::Set{Boundary},
    bc::AbstractVector{Float64};
    qdim::Int64 = 1
)
    return assemble_dirichletcondition_rhs!(
        A,
        rhs,
        extract_nodes(DI),
        bc,
        qdim = qdim
    )
end

"""
    assemble_dirichletprojection(
        ncoeffs::Int64,
        DI::Set{Int64};
        qdim::Int64 = 1
    ) -> SparseMatrixCSC{Float64, Int64}

Returns the projection matrix onto the Dirichlet nodes 
where the input ncoeffs is understood as qdim*nnodes.
"""
function assemble_dirichletprojection(
    ncoeffs::Int64,
    DI::Set{Int64};
    qdim::Int64 = 1
)
    bcids = Set{Int64}()
    for r = 1:qdim
        union!(bcids, qdim * (collect(DI) .- 1) .+ r)
    end
    nbcids = qdim*length(DI)
    
    return sparse(1:nbcids, collect(bcids), ones(nbcids), nbcids, ncoeffs)
end
