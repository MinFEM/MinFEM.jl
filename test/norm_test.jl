using MinFEM

function test_norms()
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

    mesh = import_mesh("test_square_v4.msh")
    fh = evaluate_mesh_function(mesh, f)
    abs(pnorm(1.0, fh, mesh, order=1) - 1/2) >  1e-15 && return false
    abs(pnorm(2.0, fh, mesh, order=2) - 3^(-1/2)) >  1e-15 && return false
    abs(pnorm(3.0, fh, mesh, order=3) - 4^(-1/3)) >  1e-15 && return false
    abs(pnorm(4.0, fh, mesh, order=4) - 5^(-1/4)) >  1e-15 && return false
    abs(pnorm(5.0, fh, mesh, order=5) - 6^(-1/5)) >  1e-15 && return false
    abs(pnorm(6.0, fh, mesh, order=6) - 7^(-1/6)) >  1e-15 && return false
    abs(pnorm(7.0, fh, mesh, order=7) - 8^(-1/7)) >  1e-15 && return false
    abs(pnorm(8.0, fh, mesh, order=8) - 9^(-1/8)) >  1e-15 && return false
    abs(pnorm(Inf, fh, mesh, order=1) - 1) >  1e-15 && return false

    abs(qnorm(3.0, fh, mesh, order=3) - (2/5)^(2/3)) >  1e-5 && return false
    abs(qnorm(4.0, fh, mesh, order=3) - (3/7)^(3/4)) >  1e-5 && return false
    abs(qnorm(5.0, fh, mesh, order=3) - (4/9)^(4/5)) >  1e-5 && return false
    abs(qnorm(6.0, fh, mesh, order=3) - (5/11)^(5/6)) >  1e-5 && return false
    abs(qnorm(7.0, fh, mesh, order=3) - (6/13)^(6/7)) >  1e-5 && return false
    abs(qnorm(8.0, fh, mesh, order=3) - (7/15)^(7/8)) >  1e-5 && return false

    belems = extract_elements(select_boundaries(mesh))
    abs(pnorm_boundary(1.0, fh, mesh, boundaryElements=belems, order=1) - 2) >  1e-15 && return false
    abs(pnorm_boundary(2.0, fh, mesh, boundaryElements=belems, order=2) - (2/3+1)^(1/2)) >  1e-15 && return false
    abs(pnorm_boundary(3.0, fh, mesh, boundaryElements=belems, order=3) - (2/4+1)^(1/3)) >  1e-15 && return false
    abs(pnorm_boundary(4.0, fh, mesh, boundaryElements=belems, order=4) - (2/5+1)^(1/4)) >  1e-15 && return false
    abs(pnorm_boundary(5.0, fh, mesh, boundaryElements=belems, order=5) - (2/6+1)^(1/5)) >  1e-15 && return false
    abs(pnorm_boundary(6.0, fh, mesh, boundaryElements=belems, order=6) - (2/7+1)^(1/6)) >  1e-15 && return false
    abs(pnorm_boundary(7.0, fh, mesh, boundaryElements=belems, order=7) - (2/8+1)^(1/7)) >  1e-15 && return false
    abs(pnorm_boundary(8.0, fh, mesh, boundaryElements=belems, order=8) - (2/9+1)^(1/8)) >  1e-15 && return false
    abs(pnorm_boundary(9.0, fh, mesh, boundaryElements=belems, order=9) - (2/10+1)^(1/9)) >  1e-15 && return false

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

@test test_norms()