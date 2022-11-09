# Public Documentation

---

```@meta
CurrentModule = MinFEMDev
```

Documentation for `MinFEMDev.jl`'s public interface.

See the [Internals](internals.md) page of the library for the documentation 
of internal types and functions.

## Contents

```@contents
Pages = ["public.md"]
Depth = 3
```

## Module
```@docs
MinFEMDev
```

## Types

```@docs
Mesh
Boundary
Domain
PDESystem
```

## Functions and Methods

### Mesh Generation
```@docs
unit_interval
unit_square
import_mesh
export_mesh
deform_mesh!
deform_mesh
```

### Type Handling
```@docs
select_boundaries
extract_elements
extract_nodes
```

### Function Discretization
```@docs
evaluate_mesh_function
```

### Mesh (Element) Properties
```@docs
jacobian
jacobian_boundary
elementvolume
elementvolume_boundary
elementdiameter
elementdiameter_boundary
outernormalvector
stripwidth
boundingbox
volume
barycenter
```

### Local Quadrature
```@docs
quadrature_points
quadrature_weights
quadrature_order
integral_over_reference_element
```

### Vector-Valued Coefficients Handling
```@docs
prolong_multivector
restrict_multivector
```

### Weight Computation
```@docs
assemble_weightmultivector
assemble_weightmultivector_boundary
```

### Finite Element Basis Functions
```@docs
phi
grad_phi
```

### FE Operator Assembly
```@docs
assemble_derivativematrix
assemble_laplacian
assemble_derivativematrix_boundary
assemble_normalderivativematrix
assemble_basismatrix
assemble_massmatrix
assemble_basismatrix_boundary
assemble_massmatrix_boundary
assemble_cubicterm
assemble_cubicderivativematrix
assemble_cubicsecondderivativematrix
assemble_elasticity
```

### Boundary Condition Handling
```@docs
assemble_dirichletcondition!
assemble_dirichletprojection
```

### FE Extensions
```@docs
pnorm
qnorm
pnorm_boundary
qnorm_boundary
conjugated_exponent
```

### PDE System Handling
```@docs
solve!
assemble!
refresh!
```

### Output
```@docs
write_to_vtk
write_to_vtk_boundary
open_vtkfile
open_vtkfile_boundary
save_vtkfile
write_pointdata_vtkfile!
write_celldata_vtkfile!
write_to_txt
read_from_txt
```

## Index

```@index
Pages = ["public.md"]
Module = ["MinFEMDev"]
Order = [:type, :function]
```
