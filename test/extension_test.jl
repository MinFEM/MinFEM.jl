using MinFEM

function test_conjugation()
    abs(conjugated_exponent(Inf) - 1.0) >  1e-15 && return false
    !isinf(conjugated_exponent(1.0)) && return false
    abs(conjugated_exponent(2.0) - 2.0) >  1e-15 && return false
    abs(conjugated_exponent(4.0) - 4/3) >  1e-15 && return false
    abs(conjugated_exponent(1.5) - 3.0) >  1e-15 && return false

    return true
end

function test_norm_edgecases()
    f(x) = -x[1]

    mesh::Mesh = import_mesh("test_line_v4.msh")
    fh = evaluate_mesh_function(mesh, f)

    failed::Bool = false

    try
        abs(pnorm(0.0, fh, mesh, order=1)) >  1e-15
    catch e
        if isa(e, DomainError) && occursin("conjugated Hölder exponent", e.msg)
            failed = true
        end
    end
    !failed && return false

    try
        abs(pnorm(2.0, fh[1:mesh.nelems-1], mesh, order=1)) >  1e-15
    catch e
        if isa(e, DomainError) && occursin("valid length", e.msg)
            failed = true
        end
    end
    !failed && return false

    belems = extract_elements(select_boundaries(mesh))
    try
        abs(pnorm_boundary(0.0, fh, mesh, boundaryElements=belems, order=1)) >  1e-15
    catch e
        if isa(e, DomainError) && occursin("conjugated Hölder exponent", e.msg)
            failed = true
        end
    end
    !failed && return false

    try
        abs(pnorm_boundary(2.0, fh[1:mesh.nboundelems-1], mesh, boundaryElements=belems, order=1)) >  1e-15
    catch e
        if isa(e, DomainError) && occursin("valid length", e.msg)
            failed = true
        end
    end
    !failed && return false

    abs(
        pnorm_boundary(2.0, fh, mesh, boundaryElements=belems, order=1) 
        - pnorm_boundary(2.0, fh, mesh, order=1)
    ) >  1e-15 && return false

    return true
end

function test_norms()
    f(x) = -x[1]
    g(x) = [-x[1], -x[1]]

    mesh::Mesh = import_mesh("test_line_v4.msh")
    fh = evaluate_mesh_function(mesh, f)
    fe = evaluate_quadrature_function(mesh, f, order=5)

    sol1(p) = (1/(p+1))^(1/p)

    abs(pnorm(1.0, fh, mesh, order=1) - sol1(1.0)) >  1e-15 && return false
    abs(pnorm(2.0, fh, mesh, order=2) - sol1(2.0)) >  1e-15 && return false
    abs(pnorm(3.0, fh, mesh, order=3) - sol1(3.0)) >  1e-15 && return false
    abs(pnorm(4.0, fh, mesh, order=4) - sol1(4.0)) >  1e-15 && return false
    abs(pnorm(5.0, fh, mesh, order=5) - sol1(5.0)) >  1e-15 && return false
    abs(pnorm(6.0, fh, mesh, order=6) - sol1(6.0)) >  1e-15 && return false
    abs(pnorm(7.0, fh, mesh, order=7) - sol1(7.0)) >  1e-15 && return false
    abs(pnorm(8.0, fh, mesh, order=8) - sol1(8.0)) >  1e-15 && return false
    abs(pnorm(9.0, fh, mesh, order=9) - sol1(9.0)) >  1e-15 && return false

    
    abs(pnorm(Inf, fh, mesh, order=1) - 1) >  1e-15 && return false

    abs(pnorm(1.0, fe, mesh) - sol1(1.0)) >  1e-15 && return false
    abs(pnorm(2.0, fe, mesh) - sol1(2.0)) >  1e-15 && return false
    abs(pnorm(3.0, fe, mesh) - sol1(3.0)) >  1e-15 && return false
    abs(pnorm(4.0, fe, mesh) - sol1(4.0)) >  1e-15 && return false
    abs(pnorm(5.0, fe, mesh) - sol1(5.0)) >  1e-15 && return false
    
    abs(pnorm(Inf, fe, mesh) - 0.996) >  1e-3 && return false

    abs(twonorm(fh, mesh, order=2) - sol1(2.0)) >  1e-15 && return false


    mesh = import_mesh("test_square_v4.msh")
    fh = evaluate_mesh_function(mesh, f)
    gh = evaluate_mesh_function(mesh, g, qdim=2)

    sol2(p) = (1/(p+1))^(1/p)

    abs(pnorm(1.0, fh, mesh, order=1) - sol2(1.0)) >  1e-15 && return false
    abs(pnorm(2.0, fh, mesh, order=2) - sol2(2.0)) >  1e-15 && return false
    abs(pnorm(3.0, fh, mesh, order=3) - sol2(3.0)) >  1e-15 && return false
    abs(pnorm(4.0, fh, mesh, order=4) - sol2(4.0)) >  1e-15 && return false
    abs(pnorm(5.0, fh, mesh, order=5) - sol2(5.0)) >  1e-15 && return false
    abs(pnorm(6.0, fh, mesh, order=6) - sol2(6.0)) >  1e-15 && return false
    abs(pnorm(7.0, fh, mesh, order=7) - sol2(7.0)) >  1e-15 && return false
    abs(pnorm(8.0, fh, mesh, order=8) - sol2(8.0)) >  1e-15 && return false
    
    abs(pnorm(Inf, fh, mesh, order=1) - 1) >  1e-15 && return false

    abs(twonorm(fh, mesh, order=2) - sol2(2.0)) >  1e-15 && return false

    sol3(p) = (2/(p+1))^(1/p)

    abs(pnorm(1.0, gh, mesh, order=1, qdim=2) - sol3(1.0)) >  1e-15 && return false
    abs(pnorm(2.0, gh, mesh, order=2, qdim=2) - sol3(2.0)) >  1e-15 && return false
    abs(pnorm(3.0, gh, mesh, order=3, qdim=2) - sol3(3.0)) >  1e-15 && return false
    abs(pnorm(4.0, gh, mesh, order=4, qdim=2) - sol3(4.0)) >  1e-15 && return false
    
    abs(pnorm(Inf, gh, mesh, order=1, qdim=2) - 1) >  1e-15 && return false

    abs(twonorm(gh, mesh, order=2, qdim=2) - sol3(2.0)) >  1e-15 && return false

    q(p) = p / (p-1)
    sol4(p) = (1/(q(p)+1))^(1/q(p))

    abs(qnorm(3.0, fh, mesh, order=3) - sol4(3.0)) >  1e-5 && return false
    abs(qnorm(4.0, fh, mesh, order=4) - sol4(4.0)) >  1e-5 && return false
    abs(qnorm(5.0, fh, mesh, order=5) - sol4(5.0)) >  1e-5 && return false
    abs(qnorm(6.0, fh, mesh, order=6) - sol4(6.0)) >  1e-5 && return false
    abs(qnorm(7.0, fh, mesh, order=7) - sol4(7.0)) >  1e-5 && return false
    abs(qnorm(8.0, fh, mesh, order=8) - sol4(8.0)) >  1e-5 && return false

    mesh = import_mesh("test_cube_v4.msh")
    fh = evaluate_mesh_function(mesh, f)

    sol5(p) = (1/(p+1))^(1/p)

    abs(pnorm(1.0, fh, mesh, order=1) - sol5(1.0)) >  1e-15 && return false
    abs(pnorm(2.0, fh, mesh, order=2) - sol5(2.0)) >  1e-15 && return false
    abs(pnorm(3.0, fh, mesh, order=3) - sol5(3.0)) >  1e-15 && return false
    abs(pnorm(4.0, fh, mesh, order=4) - sol5(4.0)) >  1e-15 && return false
    abs(pnorm(5.0, fh, mesh, order=5) - sol5(5.0)) >  1e-15 && return false
    abs(pnorm(6.0, fh, mesh, order=6) - sol5(6.0)) >  1e-15 && return false
    abs(pnorm(7.0, fh, mesh, order=7) - sol5(7.0)) >  1e-15 && return false
    
    abs(pnorm(Inf, fh, mesh, order=1) - 1) >  1e-15 && return false

    abs(twonorm(fh, mesh, order=2) - sol5(2.0)) >  1e-15 && return false

    return true
end

function test_boundary_norms()
    f(x) = -x[1]
    g(x) = [-x[1], -x[1]]

    mesh = import_mesh("test_square_v4.msh")
    fh = evaluate_mesh_function(mesh, f)
    fe = evaluate_quadrature_function_boundary(mesh, f, order=5)
    gh = evaluate_mesh_function(mesh, g, qdim=2)
    belems = extract_elements(select_boundaries(mesh))

    sol1(p) = (2/(p+1) + 1)^(1/p)

    abs(pnorm_boundary(1.0, fh, mesh, boundaryElements=belems, order=1) - sol1(1.0)) >  1e-15 && return false
    abs(pnorm_boundary(2.0, fh, mesh, boundaryElements=belems, order=2) - sol1(2.0)) >  1e-15 && return false
    abs(pnorm_boundary(3.0, fh, mesh, boundaryElements=belems, order=3) - sol1(3.0)) >  1e-15 && return false
    abs(pnorm_boundary(4.0, fh, mesh, boundaryElements=belems, order=4) - sol1(4.0)) >  1e-15 && return false
    abs(pnorm_boundary(5.0, fh, mesh, boundaryElements=belems, order=5) - sol1(5.0)) >  1e-15 && return false
    abs(pnorm_boundary(6.0, fh, mesh, boundaryElements=belems, order=6) - sol1(6.0)) >  1e-15 && return false
    abs(pnorm_boundary(7.0, fh, mesh, boundaryElements=belems, order=7) - sol1(7.0)) >  1e-15 && return false
    abs(pnorm_boundary(8.0, fh, mesh, boundaryElements=belems, order=8) - sol1(8.0)) >  1e-15 && return false
    abs(pnorm_boundary(9.0, fh, mesh, boundaryElements=belems, order=9) - sol1(9.0)) >  1e-15 && return false
    
    abs(pnorm_boundary(Inf, fh, mesh, boundaryElements=belems, order=1) - 1) >  1e-15 && return false

    abs(pnorm_boundary(1.0, fe, mesh, boundaryElements=belems) - sol1(1.0)) >  1e-15 && return false
    abs(pnorm_boundary(2.0, fe, mesh, boundaryElements=belems) - sol1(2.0)) >  1e-15 && return false
    abs(pnorm_boundary(3.0, fe, mesh, boundaryElements=belems) - sol1(3.0)) >  1e-15 && return false
    abs(pnorm_boundary(4.0, fe, mesh, boundaryElements=belems) - sol1(4.0)) >  1e-15 && return false
    abs(pnorm_boundary(5.0, fe, mesh, boundaryElements=belems) - sol1(5.0)) >  1e-15 && return false
    
    abs(pnorm_boundary(Inf, fe, mesh, boundaryElements=belems) - 1.0) >  1e-15 && return false

    abs(twonorm_boundary(fh, mesh, boundaryElements=belems, order=2) - sol1(2.0)) >  1e-15 && return false

    sol2(p) = (4/(p+1) + 2)^(1/p)

    abs(pnorm_boundary(1.0, gh, mesh, boundaryElements=belems, order=1, qdim=2) - sol2(1.0)) >  1e-15 && return false
    abs(pnorm_boundary(2.0, gh, mesh, boundaryElements=belems, order=2, qdim=2) - sol2(2.0)) >  1e-15 && return false
    abs(pnorm_boundary(3.0, gh, mesh, boundaryElements=belems, order=3, qdim=2) - sol2(3.0)) >  1e-15 && return false
    abs(pnorm_boundary(4.0, gh, mesh, boundaryElements=belems, order=4, qdim=2) - sol2(4.0)) >  1e-15 && return false
    abs(pnorm_boundary(Inf, gh, mesh, boundaryElements=belems, order=1, qdim=2) - 1) >  1e-15 && return false

    abs(twonorm_boundary(gh, mesh, boundaryElements=belems, order=2, qdim=2) - sol2(2.0)) >  1e-15 && return false

    q(p) = p / (p-1)
    sol3(p) = (2/(q(p)+1) + 1)^(1/q(p))

    abs(qnorm_boundary(3.0, fh, mesh, boundaryElements=belems, order=3) - sol3(3.0)) >  1e-5 && return false
    abs(qnorm_boundary(4.0, fh, mesh, boundaryElements=belems, order=4) - sol3(4.0)) >  1e-5 && return false
    abs(qnorm_boundary(5.0, fh, mesh, boundaryElements=belems, order=5) - sol3(5.0)) >  1e-5 && return false
    abs(qnorm_boundary(6.0, fh, mesh, boundaryElements=belems, order=6) - sol3(6.0)) >  1e-5 && return false
    abs(qnorm_boundary(7.0, fh, mesh, boundaryElements=belems, order=7) - sol3(7.0)) >  1e-5 && return false
    abs(qnorm_boundary(8.0, fh, mesh, boundaryElements=belems, order=8) - sol3(8.0)) >  1e-5 && return false

    mesh = import_mesh("test_cube_v4.msh")
    fh = evaluate_mesh_function(mesh, f)
    belems = extract_elements(select_boundaries(mesh))

    sol4(p) = (4/(p+1) + 1)^(1/p)

    abs(pnorm_boundary(1.0, fh, mesh, boundaryElements=belems, order=1) - sol4(1.0)) >  1e-15 && return false
    abs(pnorm_boundary(2.0, fh, mesh, boundaryElements=belems, order=2) - sol4(2.0)) >  1e-15 && return false
    abs(pnorm_boundary(3.0, fh, mesh, boundaryElements=belems, order=3) - sol4(3.0)) >  1e-15 && return false
    abs(pnorm_boundary(4.0, fh, mesh, boundaryElements=belems, order=4) - sol4(4.0)) >  1e-15 && return false
    abs(pnorm_boundary(5.0, fh, mesh, boundaryElements=belems, order=5) - sol4(5.0)) >  1e-15 && return false
    abs(pnorm_boundary(6.0, fh, mesh, boundaryElements=belems, order=6) - sol4(6.0)) >  1e-15 && return false
    abs(pnorm_boundary(7.0, fh, mesh, boundaryElements=belems, order=7) - sol4(7.0)) >  1e-15 && return false
    abs(pnorm_boundary(8.0, fh, mesh, boundaryElements=belems, order=8) - sol4(8.0)) >  1e-15 && return false

    abs(pnorm_boundary(Inf, fh, mesh, boundaryElements=belems, order=1) - 1) >  1e-15 && return false

    return true
end

@test test_conjugation()
@test test_norm_edgecases()
@test test_norms()
@test test_boundary_norms()
