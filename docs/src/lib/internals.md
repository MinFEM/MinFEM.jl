# Internal Documentation

---

```@meta
CurrentModule = MinFEM
```

Documentation for `MinFEM.jl`'s internal types and functions.

See the [Public](public.md) page of the library for the documentation of the public interface.

## Contents

```@contents
Pages = ["internals.md"]
Depth = 2
```

## Types

```@docs
Entity
```

## Quadrature

```@docs
gausslegendre_points
gausslegendre_weights
compute_coordinates_line
compute_coordinates_triangle
compute_coordinates_tetrahedron
compute_weights_line
compute_weights_triangle
compute_weights_tetrahedron
parentcoordinates
```

## Mesh Import

```@docs
import_mesh1
import_mesh2
import_mesh4
getDimFromGMSHElementType
getGMSHElementTypeFromDim
getParentBoundary
```

## Mesh Utility

```@docs
base_jacobian
circleratio
```

## FEM Utility

```@docs
stressTensor
strainTensor
```

## Index

```@index
Pages = ["internals.md"]
Module = ["MinFEM"]
Order = [:type, :function]
```
