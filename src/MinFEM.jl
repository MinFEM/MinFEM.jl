module MinFEM

using WriteVTK, SparseArrays, LinearAlgebra

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
       unit_square,
       getCellVolumes,
       open_vtk_file,
       save_vtk_file,
       write_point_data,
       write_cell_data

include("fe_types.jl")
include("mesh_methods.jl")
include("fe_methods.jl")

end # module
