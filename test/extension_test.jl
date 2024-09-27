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

    abs(pnorm(1.0, fh, mesh, order=1) - 1/2) >  1e-15 && return false
    abs(pnorm(2.0, fh, mesh, order=2) - 3^(-1/2)) >  1e-15 && return false
    abs(pnorm(3.0, fh, mesh, order=3) - 4^(-1/3)) >  1e-15 && return false
    abs(pnorm(4.0, fh, mesh, order=4) - 5^(-1/4)) >  1e-15 && return false
    abs(pnorm(5.0, fh, mesh, order=5) - 6^(-1/5)) >  1e-15 && return false
    abs(pnorm(6.0, fh, mesh, order=6) - 7^(-1/6)) >  1e-15 && return false
    abs(pnorm(7.0, fh, mesh, order=7) - 8^(-1/7)) >  1e-15 && return false
    abs(pnorm(8.0, fh, mesh, order=8) - 9^(-1/8)) >  1e-15 && return false
    abs(pnorm(9.0, fh, mesh, order=9) - 10^(-1/9)) >  1e-15 && return false
    abs(pnorm(Inf, fh, mesh, order=1) - 1) >  1e-15 && return false

    abs(pnorm(1.0, fe, mesh) - 1/2) >  1e-15 && return false
    abs(pnorm(2.0, fe, mesh) - 3^(-1/2)) >  1e-15 && return false
    abs(pnorm(3.0, fe, mesh) - 4^(-1/3)) >  1e-15 && return false
    abs(pnorm(4.0, fe, mesh) - 5^(-1/4)) >  1e-15 && return false
    abs(pnorm(5.0, fe, mesh) - 6^(-1/5)) >  1e-15 && return false
    abs(pnorm(Inf, fe, mesh) - 0.996) >  1e-3 && return false

    abs(twonorm(fh, mesh, order=2) - 3^(-1/2)) >  1e-15 && return false

    mesh = import_mesh("test_square_v4.msh")
    fh = evaluate_mesh_function(mesh, f)
    gh = evaluate_mesh_function(mesh, g, qdim=2)

    abs(pnorm(1.0, fh, mesh, order=1) - 1/2) >  1e-15 && return false
    abs(pnorm(2.0, fh, mesh, order=2) - 3^(-1/2)) >  1e-15 && return false
    abs(pnorm(3.0, fh, mesh, order=3) - 4^(-1/3)) >  1e-15 && return false
    abs(pnorm(4.0, fh, mesh, order=4) - 5^(-1/4)) >  1e-15 && return false
    abs(pnorm(5.0, fh, mesh, order=5) - 6^(-1/5)) >  1e-15 && return false
    abs(pnorm(6.0, fh, mesh, order=6) - 7^(-1/6)) >  1e-15 && return false
    abs(pnorm(7.0, fh, mesh, order=7) - 8^(-1/7)) >  1e-15 && return false
    abs(pnorm(8.0, fh, mesh, order=8) - 9^(-1/8)) >  1e-15 && return false
    abs(pnorm(Inf, fh, mesh, order=1) - 1) >  1e-15 && return false

    abs(pnorm(1.0, gh, mesh, order=1, qdim=2) - 2^(-1/2)) >  1e-15 && return false
    abs(pnorm(2.0, gh, mesh, order=2, qdim=2) - (2/3)^(1/2)) >  1e-15 && return false
    abs(pnorm(3.0, gh, mesh, order=3, qdim=2) - 2^(-1/6)) >  1e-15 && return false
    abs(pnorm(4.0, gh, mesh, order=4, qdim=2) - (4/5)^(1/4)) >  1e-15 && return false
    abs(pnorm(Inf, gh, mesh, order=1, qdim=2) - 2^(1/2)) >  1e-15 && return false

    abs(qnorm(3.0, fh, mesh, order=3) - (2/5)^(2/3)) >  1e-5 && return false
    abs(qnorm(4.0, fh, mesh, order=4) - (3/7)^(3/4)) >  1e-5 && return false
    abs(qnorm(5.0, fh, mesh, order=5) - (4/9)^(4/5)) >  1e-5 && return false
    abs(qnorm(6.0, fh, mesh, order=6) - (5/11)^(5/6)) >  1e-5 && return false
    abs(qnorm(7.0, fh, mesh, order=7) - (6/13)^(6/7)) >  1e-5 && return false
    abs(qnorm(8.0, fh, mesh, order=8) - (7/15)^(7/8)) >  1e-5 && return false

    abs(twonorm(fh, mesh, order=2) - 3^(-1/2)) >  1e-15 && return false
    abs(twonorm(gh, mesh, order=2, qdim=2) - (2/3)^(1/2)) >  1e-15 && return false

    mesh = import_mesh("test_cube_v4.msh")
    fh = evaluate_mesh_function(mesh, f)

    abs(pnorm(1.0, fh, mesh, order=1) - 1/2) >  1e-15 && return false
    abs(pnorm(2.0, fh, mesh, order=2) - 3^(-1/2)) >  1e-15 && return false
    abs(pnorm(3.0, fh, mesh, order=3) - 4^(-1/3)) >  1e-15 && return false
    abs(pnorm(4.0, fh, mesh, order=4) - 5^(-1/4)) >  1e-15 && return false
    abs(pnorm(5.0, fh, mesh, order=5) - 6^(-1/5)) >  1e-15 && return false
    abs(pnorm(6.0, fh, mesh, order=6) - 7^(-1/6)) >  1e-15 && return false
    abs(pnorm(7.0, fh, mesh, order=7) - 8^(-1/7)) >  1e-15 && return false
    abs(pnorm(Inf, fh, mesh, order=1) - 1) >  1e-15 && return false

    abs(twonorm(fh, mesh, order=2) - 3^(-1/2)) >  1e-15 && return false

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

    abs(pnorm_boundary(1.0, fh, mesh, boundaryElements=belems, order=1) - 2) >  1e-15 && return false
    abs(pnorm_boundary(2.0, fh, mesh, boundaryElements=belems, order=2) - (2/3+1)^(1/2)) >  1e-15 && return false
    abs(twonorm_boundary(fh, mesh, boundaryElements=belems, order=2) - (2/3+1)^(1/2)) >  1e-15 && return false
    abs(pnorm_boundary(3.0, fh, mesh, boundaryElements=belems, order=3) - (2/4+1)^(1/3)) >  1e-15 && return false
    abs(pnorm_boundary(4.0, fh, mesh, boundaryElements=belems, order=4) - (2/5+1)^(1/4)) >  1e-15 && return false
    abs(pnorm_boundary(5.0, fh, mesh, boundaryElements=belems, order=5) - (2/6+1)^(1/5)) >  1e-15 && return false
    abs(pnorm_boundary(6.0, fh, mesh, boundaryElements=belems, order=6) - (2/7+1)^(1/6)) >  1e-15 && return false
    abs(pnorm_boundary(7.0, fh, mesh, boundaryElements=belems, order=7) - (2/8+1)^(1/7)) >  1e-15 && return false
    abs(pnorm_boundary(8.0, fh, mesh, boundaryElements=belems, order=8) - (2/9+1)^(1/8)) >  1e-15 && return false
    abs(pnorm_boundary(9.0, fh, mesh, boundaryElements=belems, order=9) - (2/10+1)^(1/9)) >  1e-15 && return false
    abs(pnorm_boundary(Inf, fh, mesh, boundaryElements=belems, order=1) - 1) >  1e-15 && return false

    abs(pnorm_boundary(1.0, fe, mesh, boundaryElements=belems) - 2) >  1e-15 && return false
    abs(pnorm_boundary(2.0, fe, mesh, boundaryElements=belems) - (2/3+1)^(1/2)) >  1e-15 && return false
    abs(pnorm_boundary(3.0, fe, mesh, boundaryElements=belems) - (2/4+1)^(1/3)) >  1e-15 && return false
    abs(pnorm_boundary(4.0, fe, mesh, boundaryElements=belems) - (2/5+1)^(1/4)) >  1e-15 && return false
    abs(pnorm_boundary(5.0, fe, mesh, boundaryElements=belems) - (2/6+1)^(1/5)) >  1e-15 && return false
    abs(pnorm_boundary(Inf, fe, mesh, boundaryElements=belems) - 1.0) >  1e-15 && return false

    abs(pnorm_boundary(1.0, gh, mesh, boundaryElements=belems, order=1, qdim=2) - 2*(2^(1/2))) >  1e-15 && return false
    abs(pnorm_boundary(2.0, gh, mesh, boundaryElements=belems, order=2, qdim=2) - (4/3+2)^(1/2)) >  1e-15 && return false
    abs(twonorm_boundary(gh, mesh, boundaryElements=belems, order=2, qdim=2) - (4/3+2)^(1/2)) >  1e-15 && return false
    abs(pnorm_boundary(3.0, gh, mesh, boundaryElements=belems, order=3, qdim=2) - (2^(1/2)+2^(3/2))^(1/3)) >  1e-15 && return false
    abs(pnorm_boundary(4.0, gh, mesh, boundaryElements=belems, order=4, qdim=2) - (8/5+4)^(1/4)) >  1e-15 && return false
    abs(pnorm_boundary(Inf, gh, mesh, boundaryElements=belems, order=1, qdim=2) - 2^(1/2)) >  1e-15 && return false

    abs(qnorm_boundary(3.0, fh, mesh, boundaryElements=belems, order=3) - (4/5+1)^(1-(1/3))) >  1e-5 && return false
    abs(qnorm_boundary(4.0, fh, mesh, boundaryElements=belems, order=4) - (6/7+1)^(1-(1/4))) >  1e-5 && return false
    abs(qnorm_boundary(5.0, fh, mesh, boundaryElements=belems, order=5) - (8/9+1)^(1-(1/5))) >  1e-5 && return false
    abs(qnorm_boundary(6.0, fh, mesh, boundaryElements=belems, order=6) - (10/11+1)^(1-(1/6))) >  1e-5 && return false
    abs(qnorm_boundary(7.0, fh, mesh, boundaryElements=belems, order=7) - (12/13+1)^(1-(1/7))) >  1e-5 && return false
    abs(qnorm_boundary(8.0, fh, mesh, boundaryElements=belems, order=8) - (14/15+1)^(1-(1/8))) >  1e-5 && return false

    mesh = import_mesh("test_cube_v4.msh")
    fh = evaluate_mesh_function(mesh, f)
    belems = extract_elements(select_boundaries(mesh))

    abs(pnorm_boundary(1.0, fh, mesh, boundaryElements=belems, order=1) - (4/2+1)^(1/1)) >  1e-15 && return false
    abs(pnorm_boundary(2.0, fh, mesh, boundaryElements=belems, order=2) - (4/3+1)^(1/2)) >  1e-15 && return false
    abs(pnorm_boundary(3.0, fh, mesh, boundaryElements=belems, order=3) - (4/4+1)^(1/3)) >  1e-15 && return false
    abs(pnorm_boundary(4.0, fh, mesh, boundaryElements=belems, order=4) - (4/5+1)^(1/4)) >  1e-15 && return false
    abs(pnorm_boundary(5.0, fh, mesh, boundaryElements=belems, order=5) - (4/6+1)^(1/5)) >  1e-15 && return false
    abs(pnorm_boundary(6.0, fh, mesh, boundaryElements=belems, order=6) - (4/7+1)^(1/6)) >  1e-15 && return false
    abs(pnorm_boundary(7.0, fh, mesh, boundaryElements=belems, order=7) - (4/8+1)^(1/7)) >  1e-15 && return false
    abs(pnorm_boundary(8.0, fh, mesh, boundaryElements=belems, order=8) - (4/9+1)^(1/8)) >  1e-15 && return false
    abs(pnorm_boundary(Inf, fh, mesh, boundaryElements=belems, order=1) - 1) >  1e-15 && return false

    return true
end

@test test_conjugation()
@test test_norm_edgecases()
@test test_norms()
@test test_boundary_norms()
