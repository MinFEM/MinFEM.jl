"""
    MinFEM

A minimal finite element tool for demonstration and teaching in julia.

This package imports the following packages:
$(IMPORTS)
"""
module MinFEM

    using LinearAlgebra, SparseArrays, WriteVTK
    using DocStringExtensions

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
            update_mesh!,
            deform_mesh!,
            deform_mesh

    export  select_boundaries,
            select_domains,
            extract_elements,
            extract_nodes

    export  evaluate_mesh_function,
            evaluate_function

    export  jacobian,
            jacobian_boundary,
            elementvolume,
            elementvolume_boundary,
            elementbarycenter,
            elementdiameter,
            elementdiameter_boundary,
            elementratio,
            elementangle,
            outernormalvector,
            stripwidth,
            boundingbox,
            circumscribedball,
            inscribedball,
            gridsize,
            shaperegularity,
            quasiuniformity,
            volume,
            barycenter

    export  quadrature_points,
            quadrature_points_boundary,
            quadrature_weights,
            quadrature_weights_boundary,
            quadrature_order,
            parentcoordinates,
            integral_over_reference_element
    
    export  evaluate_quadrature_function,
            evaluate_quadrature_function_boundary

    export  prolong_multivector,
            restrict_multivector,
            norm_multivector

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
            assemble_dirichletcondition_rhs!,
            assemble_dirichletprojection

    export  pnorm,
            qnorm,
            pnorm_boundary,
            qnorm_boundary,
            conjugated_exponent,
            twonorm,
            twonorm_boundary

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
