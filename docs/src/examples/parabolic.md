# A Parabolic Problem

Consider the heat equation with constant heating ``f=1`` as an example for a time
dependent PDE.

```math
\begin{aligned}
u_t -\Delta u = f &\quad \text{in}\; \Omega \times (0,T)\\
u = 0 &\quad \text{on}\; \partial\Omega \times (0,T)\\
u = 0 &\quad \text{on}\; \Omega \times \{ 0\}.
\end{aligned}
```

We start by discretizing the problem uniformly in time with a finite difference approach.
Here, we include ``\theta \in [0,1]`` in order to mix implicit and explicit time stepping,
where ``\theta=1`` means fully implicit and ``\theta=0`` fully explicit.

```math
\frac{u^{k}-u^{k-1}}{\Delta t} - \Delta(\theta u^{k} + (1-\theta )u^{k-1}) = f
```

Then, for every time step ``k \in \mathbb{N}``, we need to solve the variational problem: 
Find ``u^{k} \in H_0^1(\Omega)`` such that
```math
\int_\Omega \frac{u^{k}-u^{k-1}}{\Delta t} v\, \mathrm{d}x +
\int_\Omega \nabla (\theta u^{k} + (1-\theta )u^{k-1}) \cdot \nabla v\, \mathrm{d}x =
\int_\Omega f v\, \mathrm{d}x
```
for all ``v \in H_0^1(\Omega)``.

We now multiply by ``\Delta t`` and move every known quantity to the right side
```math
\int_\Omega u^{k} v\, \mathrm{d}x + 
\theta\Delta t \int_\Omega \nabla u^{k} \cdot \nabla v\, \mathrm{d}x =
\int_\Omega u^{k-1} v\, \mathrm{d}x -
(1-\theta)\Delta t \int_\Omega \nabla u^{k-1} \cdot \nabla v\, \mathrm{d}x +
\int_\Omega f v\, \mathrm{d}x.
```
leading to the discrete formulation in terms of finite elements reading
```math
(M + \theta\Delta t L) u^{k} = (M+(\theta-1)\Delta t L) u^{k-1} + Mf.
```

The following function implements a time loop to solve the problem on the interval
``[0,T]`` in `tsteps` timesteps and outputs a *.vtu*-file for each step to a folder:

```julia
using MinFEMDev

function parabolic(mesh::Mesh, boundaryIndices::Set{Int64},
                    L::AbstractMatrix, M::AbstractMatrix,
                    f::AbstractVector, u::AbstractVector,
                    T::Float64, tsteps::Int; theta=1.0)

    dt = (T-0)/tsteps

    path = "parabolic_problem/"
    mkpath(path)
    write_to_vtk(u, mesh, "u", path*"parabolic_"*lpad(string(0), 3, '0')*".vtu")

    for i = 1:tsteps
        pde = PDESystem(A=(M + theta*dt*L), b=(M - (1.0-theta)*dt*L)*u + dt*M*f,
                        bc=zeros(mesh.nnodes), DI=boundaryIndices)
        solve!(pde)

        u = copy(pde.state)

        write_to_vtk(u, mesh, "u", path*"parabolic_"*lpad(string(i), 3, '0')*".vtu")
    end

end
```

The rest of the code is then again similar to the Poisson problem:

```julia
mesh = import_mesh("../meshes/Lshaped.msh")

boundary = select_boundaries(mesh)
boundary_nodes = extract_nodes(boundary)

source(x) = 1.0
f = evaluate_mesh_function(mesh, source)

initial_condition(x) = 0.0
u = evaluate_mesh_function(mesh, initial_condition)

L = assemble_laplacian(mesh)
M = assemble_massmatrix(mesh)

parabolic(mesh, boundary_nodes, L, M, f, u, 1.5, 100, theta=1.0)
```


From Paraview, you can access the folder created and open the file parabolic\_\*..vtu,
which will automatically load all time steps.
You can then run it via the *Play* button, as also explained in the
[tutorial](../paraview.md).
The final visualization should look similar to the following:

![Result](../assets/examples/result_parabolic.gif)
