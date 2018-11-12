module MinFEM

using SparseArrays, LinearAlgebra, WriteVTK

export PDESystem, 
       solve, 
       refresh, 
       evaluateMeshFunction, 
       asmBoundarySource,
       getDirichletProjection,
       asmLaplacian, 
       asmMassMatrix, 
       asmBoundaryMassMatrix,
       asmCubicSecondDerivativeMatrix,
       asmCubicDerivativeMatrix,
       asmCubicTerm,
       asmElasticity,
       asmDirichletCondition, 
       L2norm,
       asmGradient,
       computeGradient,
       Mesh, 
       import_mesh, 
       write_vtk_mesh

include("fe_types.jl")
include("mesh_methods.jl")
include("fe_methods.jl")

end # module
