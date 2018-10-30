quad1DW = [0.2777777777777778, 0.4444444444444444, 0.2777777777777778];

quad1DX = [0.1127016653792583, 0.5, 0.8872983346207417];

#quadX = [[0.659027622374092,  0.231933368553031],
#         [0.659027622374092,  0.109039009072877],
#         [0.231933368553031,  0.659027622374092],
#         [0.231933368553031,  0.109039009072877],
#         [0.109039009072877,  0.659027622374092],
#         [0.109039009072877,  0.231933368553031]];

#quadW = 1.0/12.0*ones(6,1);

quadX = [[1.88409405952072e-01,   7.87659461760847e-01],
         [5.23979067720101e-01,   4.09466864440735e-01],
         [8.08694385677670e-01,   8.85879595127039e-02],
         [1.06170269119576e-01,   7.87659461760847e-01],
         [2.95266567779633e-01,   4.09466864440735e-01],
         [4.55706020243648e-01,   8.85879595127039e-02],
         [2.39311322870805e-02,   7.87659461760847e-01],
         [6.65540678391645e-02,   4.09466864440735e-01],
         [1.02717654809626e-01,   8.85879595127039e-02]];

quadW = [1.93963833059595e-02,
         6.36780850998851e-02,
         5.58144204830443e-02,
         3.10342132895352e-02,
         1.01884936159816e-01,
         8.93030727728709e-02,
         1.93963833059595e-02,
         6.36780850998851e-02,
         5.58144204830443e-02];

PDESystem(A,b) = PDESystem(A,b,[],Set{Int64}(),1,[],spzeros(0,0),spzeros(0,0),[],[])
PDESystem(A,b,bc,DI) = PDESystem(A,b,bc,DI,1,[],spzeros(0,0),spzeros(0,0),[],[])
PDESystem(A,b,bc,DI,qdim) = PDESystem(A,b,bc,DI,qdim,[],spzeros(0,0),spzeros(0,0),[],[])

function assembly(S::PDESystem)
  if S.Factors == []
    if S.DI != Set{Int64}()
      vec_ind = Set{Int64}()
      for d=1:S.qdim
        union!(vec_ind, S.qdim*(collect(S.DI).-1).+d)
      end
      ii = length(vec_ind)
      m,n = size(S.A)
      S.SystemMatrix = spzeros(n+ii, n+ii)
      S.B = sparse(1:ii,collect(vec_ind),ones(ii),ii,m)
      S.SystemMatrix[1:n, 1:n] = S.A
      S.SystemMatrix[1:n, (n+1):end] = S.B'
      S.SystemMatrix[(n+1):end, 1:n] = S.B
    else
      S.SystemMatrix = S.A
    end
    S.Factors = lu(S.SystemMatrix)
  end
end

function refresh(S::PDESystem)
  S.Factors = []
end

"""
    solve(S::PDESystem)

First tries to set up the system matrix with multipliers for Dirichlet conditions.
If the system has already been used before, this step is skipped.
This is determined depending on an existing factorization of the system matrix.
If the stiffness matrix or Dirichlet conditions have changes, one should invole refresh() first.

Finally the system is solved via matrix factorization.
"""
function solve(S::PDESystem)
  assembly(S)
  vec_ind = Set{Int64}()
  for d=1:S.qdim
    union!(vec_ind, S.qdim*(collect(S.DI).-1).+d)
  end
  S.rhs = [S.b;S.bc[collect(vec_ind)]]
  S.state = (S.Factors\S.rhs)[1:length(S.b)]
end

"""
    evaluateMeshFunction(mesh::Mesh, f::Function; region=[], qdim=1)

Evaluate a given function object f on all nodes of the mesh. Or, if region is specified,
it is assumed to be a set or vector of node indices, where f should be evaluated.
"""
function evaluateMeshFunction(mesh::Mesh, f::Function; region=[], qdim=1)
  v = zeros(Float64, qdim*mesh.nnodes)
  nodes = 1:mesh.nnodes
  if region != []
    nodes = collect(region)
  end
  if(qdim==1)
    for i in nodes
      v[i] = f(mesh.Nodes[i])
    end
  else
    for i in nodes
      v[qdim*(i-1)+1:qdim*i] = f(mesh.Nodes[i])
    end
  end
  return v
end

function Phi(ii::Int64, x::AbstractVector)
  if ii==1
    return 1.0-x[1]-x[2]
  elseif ii==2
    return x[1]
  elseif ii==3
    return x[2]
  end
end

function PhiEdge(ii::Int64, x::Float64)
  if ii==1
    return 1.0-x
  elseif ii==2
    return x
  end
end

function gradPhi(ii::Int64)
  grad = zeros(2,1)
  if ii==1
    grad = [-1.0; -1.0]
  elseif ii==2
    grad = [1.0; 0.0]
  elseif ii==3
    grad = [0.0; 1.0]
  end
  return grad
end

"""
    function asmDirichletCondition(SM, DI::Set{Int64}, rhs=[], bc=[]; qdim=1)

Modify a stiffness matrix and a right hand side according to the given Dirichlet conditions.
DI has to be the set of node indices for which the condition should be active.
For vector valued states either DI can be set to each component that should have a
Dirichlet condtion or qdim is set, if all components should have the condition.
"""
function asmDirichletCondition(SM, DI::Set{Int64}, rhs=[], bc=[]; qdim=1)
  if rhs != [] && bc != []
    for i in DI
      for d=1:qdim
        ii = qdim*(i-1)+d
        bcind = SM[:, ii].nzind
        rhs[bcind] -= SM[bcind, ii]*bc[ii]
      end
    end
    for i in DI
      for d=1:qdim
        ii = qdim*(i-1)+d
        rhs[ii] = bc[ii]
      end
    end
  end
  for i in DI
    for d=1:qdim
      ii = qdim*(i-1)+d
      SM[ii, SM[ii,:].nzind] *= 0.0;
      SM[:, ii] *= 0.0;
      SM[ii, ii] = 1.0;
    end
  end
end

function asmSparseMatrix(mesh::Mesh, D::Dict{Tuple{Int64,Int64}, Float64}; nrows=-1, ncols=-1)
  nn = length(D)
  II = zeros(Int64, nn)
  JJ = zeros(Int64, nn)
  AA = zeros(Float64, nn)
  for (i, pair) in enumerate(D)
    II[i] = pair[1][1]
    JJ[i] = pair[1][2]
    AA[i] = pair[2]
  end
  return sparse(II,JJ,AA, maximum([nrows,maximum(II)]), maximum([ncols,maximum(JJ)]))
end


"""
    asmLaplacian(mesh::Mesh)

Assemble the Laplacian stiffness matrix for all elements in the mesh.
"""
function asmLaplacian(mesh::Mesh)
  D = Dict{Tuple{Int64,Int64}, Float64}()

  for el=1:mesh.nelems
    nodes = mesh.Triangles[el]
    (detJ, J) = Jacobian(mesh, el)
    elemMat = zeros(3,3)
    for i=1:3
      for j=1:3
        elemMat[i,j] += (J*gradPhi(i))'*(J*gradPhi(j))
      end
    end
    elemMat *= sum(quadW)*detJ

    for i=1:3
      for j=1:3
        if(in((nodes[i], nodes[j]), keys(D)))
          D[(nodes[i], nodes[j])] += elemMat[i,j]
        else
          D[(nodes[i], nodes[j])] = elemMat[i,j]
        end
      end
    end
  end

  return asmSparseMatrix(mesh, D)
end

"""
    asmMassMatrix(mesh::Mesh; qdim=1)

Assemble a mass matrix for all elements of the given mesh.
"""
function asmMassMatrix(mesh::Mesh; qdim=1)
  D = Dict{Tuple{Int64,Int64}, Float64}()

  for el=1:mesh.nelems
    nodes = mesh.Triangles[el]
    (detJ, J) = Jacobian(mesh, el)
    elemMat = zeros(3,3)
    for i=1:3
      for j=1:3
        for (q, x) in enumerate(quadX)
          elemMat[i,j] += Phi(i, x) * Phi(j, x) * quadW[q] * detJ
        end
      end
    end

    for d=1:qdim
      for i=1:3
        for j=1:3
          ii = qdim*(nodes[i]-1)+d
          jj = qdim*(nodes[j]-1)+d
          if(in((ii, jj), keys(D)))
            D[(ii, jj)] += elemMat[i,j]
          else
            D[(ii, jj)] = elemMat[i,j]
          end
        end
      end
    end
  end

  return asmSparseMatrix(mesh, D)
end


"""
    asmBoundaryMassMatrix(mesh::Mesh, BoundaryEdges=Set{Int64}(-1); qdim=1)

Assemble a mass matrix for the given set of boundary edges.
"""
function asmBoundaryMassMatrix(mesh::Mesh, BoundaryEdges=Set{Int64}(-1); qdim=1)
  D = Dict{Tuple{Int64,Int64}, Float64}()

  if(in(-1, BoundaryEdges))
    BoundaryEdges = 1:mesh.nedges
  end

  for el in BoundaryEdges
    edge = mesh.Edges[el]
    detJ = EdgeJacobian(mesh, el)
    elemMat = zeros(2,2)
    for i=1:2
      for j=1:2
        for (q, x) in enumerate(quad1DX)
          elemMat[i,j] += PhiEdge(i, x) * PhiEdge(j, x) * quad1DW[q] * detJ
        end
      end
    end

    for d=1:qdim
      for i=1:2
        for j=1:2
          ii = qdim*(edge[i]-1)+d
          jj = qdim*(edge[j]-1)+d
          if(in((ii, jj), keys(D)))
            D[(ii, jj)] += elemMat[i,j]
          else
            D[(ii, jj)] = elemMat[i,j]
          end
        end
      end
    end
  end

  return asmSparseMatrix(mesh,D,nrows=qdim*mesh.nnodes,ncols=qdim*mesh.nnodes)
end


"""
    asmBoundarySource(mesh::Mesh, S::Array{Float64,1}, BoundaryEdges=Set{Int64}(-1); qdim=1)

Assemble the source S on the given set of boundary edges.
"""
function asmBoundarySource(mesh::Mesh, S::Array{Float64,1}, BoundaryEdges=Set{Int64}(-1); qdim=1)
  V = zeros(mesh.nnodes)

  if(in(-1, BoundaryEdges))
    BoundaryEdges = 1:mesh.nedges
  end

  for el in BoundaryEdges
    edge = mesh.Edges[el]
    detJ = EdgeJacobian(mesh, el)
    for i=1:2
      for j=1:2
        for (q, x) in enumerate(quad1DX)
          V[edge[i]] += S[edge[j]]*PhiEdge(j, x) * PhiEdge(i, x) * quad1DW[q] * detJ
        end
      end
    end
  end

  return V
end

"""
    L2norm(mesh::Mesh, Mass::SparseMatrixCSC{Float64,Int64}, v::AbstractVector; qdim=1)

Computes the L2 norm of a vector v according to the mass matrix Mass.
"""
function L2norm(Mass::SparseMatrixCSC{Float64,Int64}, v::AbstractVector; qdim=1)
  return sqrt(v'*Mass*v)
end

"""
    asmCubicTerm(mesh::Mesh, y::AbstractVector)

The cubic term of the standard semilinear parabolic equation.
"""
function asmCubicTerm(mesh::Mesh, y::AbstractVector)
  V = zeros(mesh.nnodes)
  for el=1:mesh.nelems
    nodes = mesh.Triangles[el]
    (detJ, J) = Jacobian(mesh, el)
    y_cubic = zeros(length(quadW))
    for i=1:3
      for (q, x) in enumerate(quadX)
        y_cubic[q] += (y[nodes[i]]*Phi(i, x))^3
      end
    end

    for i=1:3
      for (q, x) in enumerate(quadX)
        V[nodes[i]] += y_cubic[q] * Phi(i, x) * quadW[q] * detJ
      end
    end

  end

  return V
end

"""
    asmCubicDerivativeMatrix(mesh::Mesh, y::AbstractVector)

Assembly of the linearization of the cubic term of the standard semilinear elliptic equation.
"""
function asmCubicDerivativeMatrix(mesh::Mesh, y::AbstractVector)
  D = Dict{Tuple{Int64,Int64}, Float64}()

  for el=1:mesh.nelems
    nodes = mesh.Triangles[el]
    y_quadratic = zeros(length(quadW))
    y_quad = zeros(length(quadW))
    for i=1:3
      for (q, x) in enumerate(quadX)
        y_quad[q] += y[nodes[i]]*Phi(i, x)
      end
    end

    (detJ, J) = Jacobian(mesh, el)
    elemMat = zeros(3,3)
    for i=1:3
      for j=1:3
        for (q, x) in enumerate(quadX)
          elemMat[i,j] += 3.0*y_quad[q]^2 * Phi(i, x) * Phi(j, x) * quadW[q] * detJ
        end
      end
    end

    for i=1:3
      for j=1:3
        if(in((nodes[i], nodes[j]), keys(D)))
          D[(nodes[i], nodes[j])] += elemMat[i,j]
        else
          D[(nodes[i], nodes[j])] = elemMat[i,j]
        end
      end
    end
  end

  return asmSparseMatrix(mesh, D)
end

"""
    asmCubicSecondDerivativeMatrix(mesh::Mesh, y::AbstractVector, p::AbstractVector)

Assembly of the second derivative of the cubic term of the standard semilinear elliptic equation 
around the state (y,p).
"""
function asmCubicSecondDerivativeMatrix(mesh::Mesh, y::AbstractVector, p::AbstractVector)
  D = Dict{Tuple{Int64,Int64}, Float64}()

  for el=1:mesh.nelems
    nodes = mesh.Triangles[el]
    y_quad = zeros(length(quadW))
    p_quad = zeros(length(quadW))
    for i=1:3
      for (q, x) in enumerate(quadX)
        y_quad[q] += y[nodes[i]]*Phi(i, x)
        p_quad[q] += p[nodes[i]]*Phi(i, x)
      end
    end

    (detJ, J) = Jacobian(mesh, el)
    elemMat = zeros(3,3)
    for i=1:3
      for j=1:3
        for (q, x) in enumerate(quadX)
          elemMat[i,j] += 6.0*y_quad[q]*p_quad[q] * Phi(i, x) * Phi(j, x) * quadW[q] * detJ
        end
      end
    end

    for i=1:3
      for j=1:3
        if(in((nodes[i], nodes[j]), keys(D)))
          D[(nodes[i], nodes[j])] += elemMat[i,j]
        else
          D[(nodes[i], nodes[j])] = elemMat[i,j]
        end
      end
    end
  end

  return asmSparseMatrix(mesh, D)
end

"""
    strainTensor(grad::AbstractMatrix)

Strain tensor for linear elasticity with gradient as argument.
"""
function strainTensor(grad::AbstractMatrix)
  return 0.5*(grad + grad')
end

"""
    stressTensor(grad::AbstractMatrix, lambda::Float64, mu::Float64)

Strain tensor for linear elasticity with gradient as argument.
"""
function stressTensor(grad::AbstractMatrix, lambda::Float64, mu::Float64)
  epsilon = strainTensor(grad)
  return lambda * tr(epsilon) * [1.0 0.0;0.0 1.0] + 2.0 * mu * epsilon
end

"""
    asmElasticity(mesh::Mesh, lambda::Float64, mu::Float64)

Assembly of the linear elasticity stiffness matrix for constant coefficients lambda and mu.
"""
function asmElasticity(mesh::Mesh, lambda::Float64, mu::Float64)
  D = Dict{Tuple{Int64,Int64}, Float64}()

  qdim = 2

  for el=1:mesh.nelems
    nodes = mesh.Triangles[el]
    (detJ, J) = Jacobian(mesh, el)
    elemMat = zeros(qdim*3,qdim*3)
    for i=1:3
      for j=1:3
        for ic=1:qdim
          for jc=1:qdim
            grad_i = zeros(2,2)
            grad_i[:,ic] = J*gradPhi(i)
            grad_j = zeros(2,2)
            grad_j[:,jc] = J*gradPhi(j)
            elemMat[qdim*(i-1)+ic,qdim*(j-1)+jc] += dot(stressTensor(grad_j,lambda,mu), strainTensor(grad_i))
          end
        end
      end
    end
    elemMat *= sum(quadW)*detJ

    for i=1:3
      for j=1:3
        for ic=1:qdim
          for jc=1:qdim
            ii = qdim*(nodes[i]-1)+ic
            jj = qdim*(nodes[j]-1)+jc
            if(in((ii, jj), keys(D)))
              D[(ii, jj)] += elemMat[qdim*(i-1)+ic,qdim*(j-1)+jc]
            else
              D[(ii, jj)] = elemMat[qdim*(i-1)+ic,qdim*(j-1)+jc]
            end
          end
        end
      end
    end
  end

  return asmSparseMatrix(mesh, D)
end
