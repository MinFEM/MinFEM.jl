"""
    Boundary

Structure holding sets of node and edge indices for one particular physical boundary.
"""
mutable struct Boundary
  Nodes::Set{Int64}
  Edges::Set{Int64}
end

"""
    Mesh

Structure for a triangular finite element mesh with volume and boundary markers.
"""
mutable struct Mesh
  "number of nodes"
  nnodes::Int64

  "number of triangles"
  nelems::Int64

  "number of physical (marked) edges"
  nedges::Int64

  "list of node coordinates"
  Nodes::Array{Array{Float64,1},1}

  "list of triangle node indices"
  Triangles::Array{Array{Int64,1},1}

  "list of physical edges"
  Edges::Array{Array{Int64,1},1}

  "dictionary of marked boundaries"
  Boundaries::Dict{Int64,Boundary}

  "dictionary of marked volume regions"
  Volumes::Dict{Int64,Set{Int64}}
end

"""
  PDESystem

Structure holding all information to describe simple PDEs with Dirichlet realized
resolved by multipliers.
"""
mutable struct PDESystem
  "stiffness matrix"
  A::SparseMatrixCSC{Float64,Int64}

  "load vector"
  b::Vector{Float64}

  "Dirichlet boundary values"
  bc::Vector{Float64}

  "Dirichlet nodes"
  DI::Set{Int64}

#  "Dirichlet indices, in case of qdim!=1"
#  vec_ind::Set{Int64}

  "dimension of vector-valued state"
  qdim::Int64

  "system matrix factorization"
  Factors

  "system of stiffness matrix and Dirichlet conditinons"
  SystemMatrix::SparseMatrixCSC{Float64,Int64}

  "Dirichlet projection matrix"
  B::SparseMatrixCSC{Float64,Int64}

  "solution vector"
  state::Vector{Float64}

  "rhs vector with source term and dirichlet conditions"
  rhs::Vector{Float64}

end
