"""
    MinFEMDev

A minimal finite element tool for demonstration and teaching in julia.
"""
module MinFEMDev

using LinearAlgebra, SparseArrays, WriteVTK

include("mesh.jl")
include("quadrature.jl")
include("fem.jl")
include("extensions.jl")
include("pdesystem.jl")
include("output.jl")

export  Mesh,
        Boundary,
        Domain,
        PDESystem

export  unit_square,
        unit_interval,
        import_mesh,
        export_mesh,
        deform_mesh!,
        deform_mesh

export  select_boundaries,
        extract_elements,
        extract_nodes

export  evaluate_mesh_function

export  jacobian,
        jacobian_boundary,
        elementvolume,
        elementvolume_boundary,
        elementdiameter,
        elementdiameter_boundary,
        outernormalvector,
        stripwidth,
        boundingbox,
        volume,
        barycenter

export  quadrature_points,
        quadrature_weights,
        quadrature_order,
        integral_over_reference_element

export  prolong_multivector,
        restrict_multivector

export  assemble_weightmultivector,
        assemble_weightmultivector_boundary

export  phi,
        grad_phi

export  assemble_derivativematrix,
        assemble_laplacian,
        assemble_derivativematrix_boundary,
        assemble_normalderivativematrix,
        assemble_basismatrix,
        assemble_massmatrix,
        assemble_basismatrix_boundary,
        assemble_massmatrix_boundary,
        assemble_cubicterm,
        assemble_cubicderivativematrix,
        assemble_cubicsecondderivativematrix,
        assemble_elasticity

export  assemble_dirichletcondition!,
        assemble_dirichletprojection

export  pnorm,
        qnorm,
        pnorm_boundary,
        qnorm_boundary,
        conjugated_exponent

export  solve!,
        assemble!,
        refresh!

export  write_to_vtk,
        write_to_vtk_boundary,
        open_vtkfile,
        open_vtkfile_boundary,
        save_vtkfile,
        write_pointdata_vtkfile!,
        write_celldata_vtkfile!,
        write_to_txt,
        read_from_txt

end # module
