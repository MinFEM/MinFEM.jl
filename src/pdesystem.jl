"""
    PDESystem

Structure holding all information to describe simple PDEs 
with Dirichlet boundary conditions.
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

    "dimension of vector-valued state"
    qdim::Int64

    "system matrix factorization"
    Factors::Any

    "system of stiffness matrix and Dirichlet conditinons"
    SystemMatrix::SparseMatrixCSC{Float64,Int64}

    "Dirichlet projection matrix"
    B::SparseMatrixCSC{Float64,Int64}

    "solution vector"
    state::Vector{Float64}

    "rhs vector with source term and dirichlet conditions"
    rhs::Vector{Float64}

end

PDESystem(;A=spzeros(0,0), b=[], bc=[], DI=Set{Int64}(), qdim=1, Factors=[],
            SystemMatrix=spzeros(0,0), B=spzeros(0,0), state=[], rhs=[]) =
                PDESystem(A, b, bc, DI, qdim, Factors, SystemMatrix, B, state, rhs)

"""
    assemble!(S::PDESystem)

If the system has not been used before, sets up the system matrix 
with multipliers for Dirichlet conditions and factorizes it.
"""
function assemble!(S::PDESystem)
    if S.Factors == []
        if !isempty(S.DI)
            m,n = size(S.A)
            S.B = assemble_dirichletprojection(m, S.DI, qdim=S.qdim)
            nbcids = length(S.DI) * S.qdim

            S.SystemMatrix = [S.A                   S.B';
                              S.B spzeros(nbcids,nbcids)]
        else
            S.SystemMatrix = S.A
        end

        S.Factors = lu(S.SystemMatrix)
    end
end

"""
    refresh!(S::PDESystem)

Recomputes the factorization of the stiffness matrix using `assemble!()`.
"""
function refresh!(S::PDESystem)
    S.Factors = []
    assemble!(S)
end

"""
    solve!(S::PDESystem)

First tries to set up the system matrix with multipliers for Dirichlet conditions.
If the system has already been used before, this step is skipped.
This is determined depending on an existing factorization of the system matrix.
If the stiffness matrix or Dirichlet conditions have changes, 
one should invole refresh() first.
Finally the system is solved via matrix factorization.
"""
function solve!(S::PDESystem)
    assemble!(S)
    S.rhs = [S.b; S.B*S.bc]
    S.state = (S.Factors\S.rhs)[1:length(S.b)]
end
