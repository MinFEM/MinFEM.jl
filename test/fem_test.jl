using MinFEM, SparseArrays

function test_restriction()
    x = [1.0, 2.0 , 3.0, 4.0]

    restrict_multivector(x, 4, 1, block=1) != [1.0, 2.0, 3.0, 4.0] && return false
    restrict_multivector(x, 2, 2, block=1) != [3.0, 7.0] && return false
    restrict_multivector(x, 2, 1, block=2) != [1.0, 2.0, 3.0, 4.0] && return false
    restrict_multivector(x, 1, 2, block=2) != [3.0, 7.0] && return false

    return true
end

function test_prolongation()
    x = [1.0, 2.0]

    prolong_multivector(x, 2, 1, block=1) != [1.0, 2.0] && return false
    prolong_multivector(x, 2, 2, block=1) != [1.0, 1.0, 2.0, 2.0] && return false
    prolong_multivector(x, 2, 1, block=2) != [1.0, 2.0] && return false
    prolong_multivector(x, 1, 2, block=2) != [1.0, 2.0, 1.0, 2.0] && return false

    return true
end

function test_laplacian()
    mesh = import_mesh("test_square_v4.msh")

    L1 = assemble_laplacian(mesh, qdim=1)
    w = assemble_weightmultivector(mesh, qdim=2, order=1)
    D = assemble_derivativematrix(mesh, qdim=1)
    L2 = assemble_laplacian(D, w)
    any(abs.(nonzeros(L1) .- nonzeros(L2)) .> 1e-13) && return false

    L1 = assemble_laplacian(mesh, qdim=2)
    w = assemble_weightmultivector(mesh, qdim=2*2, order=1)
    D = assemble_derivativematrix(mesh, qdim=2)
    L2 = assemble_laplacian(D, w)
    any(abs.(nonzeros(L1) .- nonzeros(L2)) .> 1e-13) && return false

    return true
end

function test_derivative_boundary()
    mesh = import_mesh("test_square_v4.msh")
    f(x) = -x[1]
    g(x) = [-x[1],-x[1]]

    fh = evaluate_mesh_function(mesh, f)
    gh = evaluate_mesh_function(mesh, g, qdim=2)

    df(x) = [-1.0 ,0.0]
    dfh = evaluate_quadrature_function_boundary(mesh, df, qdim=2)
    dg(x) = [-1.0, -1.0, 0.0, 0.0]
    dgh = evaluate_quadrature_function_boundary(mesh, dg, qdim=4)

    D1 = assemble_derivativematrix_boundary(mesh, qdim=1)
    D2 = assemble_derivativematrix_boundary(mesh, qdim=2)

    any(abs.(D1*fh .- dfh) .> 1e-13) && return false
    any(abs.(D2*gh .- dgh) .> 1e-13) && return false

    return true
end

function test_derivative_normal()
    mesh = import_mesh("test_square_v4.msh")
    f(x) = -x[1]
    fh = evaluate_mesh_function(mesh, f)

    D = assemble_normalderivativematrix(mesh, qdim=1)

    function df(x)
        if x[1] < 1e-10
            return 1
        elseif x[1] < 1.0
            return 0
        else 
            return -1
        end
    end
    dfh = evaluate_quadrature_function_boundary(mesh, df, qdim=1)

    any(abs.(D*fh .- dfh) .> 1e-13) && return false

    return true
end

function test_mass()
    mesh = import_mesh("test_square_v4.msh")

    M1 = assemble_massmatrix(mesh, qdim=1)
    w = assemble_weightmultivector(mesh, qdim=1, order=3)
    E = assemble_basismatrix(mesh, qdim=1)
    M2 = assemble_massmatrix(E, w)
    any(abs.(nonzeros(M1) .- nonzeros(M2)) .> 1e-13) && return false

    M1 = assemble_massmatrix(mesh, qdim=2, order=5)
    w = assemble_weightmultivector(mesh, qdim=2, order=5)
    E = assemble_basismatrix(mesh, qdim=2, order=5)
    M2 = assemble_massmatrix(E, w)
    any(abs.(nonzeros(M1) .- nonzeros(M2)) .> 1e-13) && return false

    return true
end

function test_mass_boundary()
    mesh = import_mesh("test_square_v4.msh")

    M1 = assemble_massmatrix_boundary(mesh, qdim=1)
    w = assemble_weightmultivector_boundary(mesh, qdim=1, order=3)
    E = assemble_basismatrix_boundary(mesh, qdim=1)
    M2 = assemble_massmatrix_boundary(E, w)
    any(abs.(nonzeros(M1) .- nonzeros(M2)) .> 1e-13) && return false

    M1 = assemble_massmatrix_boundary(mesh, qdim=2, order=5)
    w = assemble_weightmultivector_boundary(mesh, qdim=2, order=5)
    E = assemble_basismatrix_boundary(mesh, qdim=2, order=5)
    M2 = assemble_massmatrix_boundary(E, w)
    any(abs.(nonzeros(M1) .- nonzeros(M2)) .> 1e-13) && return false

    return true
end

function test_dirichletcondition()
    mesh = import_mesh("test_square_v4.msh")

    bnd = select_boundaries(mesh, 1001, 1002)
    belemes = extract_nodes(bnd)

    L1 = assemble_laplacian(mesh)
    L2 = assemble_laplacian(mesh)
    assemble_dirichletcondition!(L1, bnd)
    assemble_dirichletcondition!(L2, belemes)

    any(abs.(nonzeros(L1) .- nonzeros(L2)) .> 1e-13) && return false
    
    return true
end

@test test_restriction()
@test test_prolongation()
@test test_laplacian()
@test test_derivative_boundary()
@test test_derivative_normal()
@test test_mass()
@test test_mass_boundary()
@test test_dirichletcondition()
