# Internal Documentation

---

```@meta
CurrentModule = MinFEM
```

Documentation for `MinFEM.jl`'s internal types and functions.

See the [Public](public.md) page of the library for the documentation of the public interface.

## Contents

```@contents
Pages = ["internal.md"]
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
gmsh_dimfromtype
gmsh_typefromdim
parentboundary
sort_boundaryelement
```

## Mesh Utility
```@docs
base_jacobian
circumscribedball1d
circumscribedball2d
circumscribedball3d
inscribedball1d
inscribedball2d
inscribedball3d
```

## FEM Utility

```@docs
stress_tensor
strain_tensor
```

## Index

```@index
Pages = ["internal.md"]
Order = [:type, :function]
```
