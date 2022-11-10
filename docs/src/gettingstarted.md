# Getting Started

MinFEM.jl is available as an registered Julia package and thus the installation 
can be done quickly using the package manager.
To do so, open the julia REPL, hit the **]** key to change into the package mode.
The line will now start with something like `(@v1.7) pkg>` where the numbers denote 
you current minor julia version, indicating that you are modifying your full installation.
In case there is a package name featured instead, 
because you managed dependencies of it before, you need to restart the full REPL.
Now you can add the package by typing

```
add MinFEM
test MinFEM
```

If you are interested in modifying the code, you might also clone the GitHub repository 
from [here](https://github.com/MinFEM/MinFEM.jl) to your local machine.
Then you can add it as a development package to your Julia installation 
by using the package tool in the REPL as before and typing

```
dev C:\your\path\to\minfem
test MinFEM
```

From now on, you can use MinFEM by including it on top of you code file

```julia
using MinFEM
```

In general the workflow starts by creating a mesh of the computational domain.
This can either be done by importing a GMSH file or using the internal functions 
for creating a unit mesh.

```julia
mesh = import_mesh("../meshes/square.msh")

mesh = unit_square(50)
```

Subsequently, you can specify physical boundary sets that are used
to specify boundary conditions later on.
The general syntax returns you a set of the type `Boundary`, which contains either
all physical boundaries specified, or just the selected indices.
The indices are either specified in the GMSH file or for the unit meshes given in
their respective documentation.
Here, 1001 and 1003 select the bottom and left boundary of the unit square.
For some applications you might need the all the nodes or elements of a boundary set. 
These can be easily extracted and are then given by an unordered set `Set{Int64}`.

```julia
boundary = select_boundaries(mesh)
boundary = select_boundaries(mesh, 1001, 1003)

boundaryNodes = extract_nodes(boundary)
boundaryNodes = extract_elements(boundary)
```

Next, operators can be assembled. 
There are multiple operators available all following the same syntax.
For a full list check the [Public Documentation](@ref) page,
here we only introduce the two most basic ones.

```julia
L = assemble_laplacian(mesh)
M = assemble_massmatrix(mesh)
```

While one could now manually code the solution of the finite element linear system,
MinFEM supports the creation of a `PDESystem`.
With this you can set up the linear system as you would in a classical sense
and automatically add boundary conditions.
Then simply solve the system and the solution will be contained in the object.

```julia
pde = PDESystem(A=L, b=M*ones(mesh.nnodes), bc=zeros(mesh.nnodes), DI=boundaryNodes)
solve!(pde)
```

Finally, you can write the solution to a file.
Most commonly, this will be a *.vtu* file to be read with Paraview for visualization.
For the storage of intermediate result it might also be useful to save it
as a simple *.txt* file.
This can be read later on and could also be used for visualization with Plots.jl.

```julia
write_to_vtk(pde.state, mesh, "u", "getting_started")

write_to_txt(pde.state, mesh, "u", "getting_started")
```

For more detailed tutorial workflows check out the Examples section of this manual.
These will include multiple operators as well as a time dependent problem,
a non-linear problem and how to handle vector-valued functions.
