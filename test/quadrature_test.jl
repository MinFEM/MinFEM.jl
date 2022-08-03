using MinFEM

function test_quadrature()
    f0(x) = 1
    f1(x) = x[1]
    f2(x) = x[1]^2
    f3(x) = x[1]^3
    f4(x) = x[1]^4
    f5(x) = x[1]^5
    f6(x) = x[1]^6
    f7(x) = x[1]^7
    f8(x) = x[1]^8
    f9(x) = x[1]^9
    f10(x) = x[1]^10

    failed::Bool = false
    try
        abs(integral_over_reference_element(f0, -1, order=0) - 1) >  1e-15
    catch e
        if isa(e, ErrorException) && occursin("Dimension -1", e.msg)
            failed = true
        end
    end
    !failed && return false

    try
        abs(integral_over_reference_element(f0, 4, order=0) - 1) >  1e-15
    catch e
        if isa(e, ErrorException) && occursin("Dimension 4", e.msg)
            failed = true
        end
    end
    !failed && return false

    failed = false
    try
        abs(integral_over_reference_element(f0, 1, order=-1) - 1) >  1e-15
    catch e
        if isa(e, ErrorException) && occursin("Order -1", e.msg)
            failed = true
        end
    end
    !failed && return false
    
    abs(integral_over_reference_element(f0, 1, order=0) - 1/1) >  1e-15 && return false
    abs(integral_over_reference_element(f1, 1, order=1) - 1/2) >  1e-15 && return false
    abs(integral_over_reference_element(f2, 1, order=2) - 1/3) >  1e-15 && return false
    abs(integral_over_reference_element(f3, 1, order=3) - 1/4) >  1e-15 && return false
    abs(integral_over_reference_element(f4, 1, order=4) - 1/5) >  1e-15 && return false
    abs(integral_over_reference_element(f5, 1, order=5) - 1/6) >  1e-15 && return false
    abs(integral_over_reference_element(f6, 1, order=6) - 1/7) >  1e-15 && return false
    abs(integral_over_reference_element(f7, 1, order=7) - 1/8) >  1e-15 && return false
    abs(integral_over_reference_element(f8, 1, order=8) - 1/9) >  1e-15 && return false
    abs(integral_over_reference_element(f9, 1, order=9) - 1/10) >  1e-15 && return false
    
    failed = false
    try
        abs(integral_over_reference_element(f10, 1, order=10) - 1/11) >  1e-15
    catch e
        if isa(e, ErrorException) && occursin("Order 10", e.msg) && occursin("1D", e.msg)
            failed = true
        end
    end
    !failed && return false

    abs(integral_over_reference_element(f0, 2, order=0) - 1/2) >  1e-15 && return false
    abs(integral_over_reference_element(f1, 2, order=1) - 1/6) >  1e-15 && return false
    abs(integral_over_reference_element(f2, 2, order=2) - 1/12) >  1e-15 && return false
    abs(integral_over_reference_element(f3, 2, order=3) - 1/20) >  1e-15 && return false
    abs(integral_over_reference_element(f4, 2, order=4) - 1/30) >  1e-15 && return false
    abs(integral_over_reference_element(f5, 2, order=5) - 1/42) >  1e-15 && return false
    abs(integral_over_reference_element(f6, 2, order=6) - 1/56) >  1e-15 && return false
    abs(integral_over_reference_element(f7, 2, order=7) - 1/72) >  1e-15 && return false
    abs(integral_over_reference_element(f8, 2, order=8) - 1/90) >  1e-15 && return false
    
    failed = false
    try
        abs(integral_over_reference_element(f9, 2, order=9) - 1/110) >  1e-15
    catch e
        if isa(e, ErrorException) && occursin("Order 9", e.msg) && occursin("2D", e.msg)
            failed = true
        end
    end
    !failed && return false

    abs(integral_over_reference_element(f0, 3, order=0) - 1/6) >  1e-15 && return false
    abs(integral_over_reference_element(f1, 3, order=1) - 1/24) >  1e-15 && return false
    abs(integral_over_reference_element(f2, 3, order=2) - 1/60) >  1e-15 && return false
    abs(integral_over_reference_element(f3, 3, order=3) - 1/120) >  1e-15 && return false
    abs(integral_over_reference_element(f4, 3, order=4) - 1/210) >  1e-15 && return false
    abs(integral_over_reference_element(f5, 3, order=5) - 1/336) >  1e-15 && return false
    abs(integral_over_reference_element(f6, 3, order=6) - 1/504) >  1e-15 && return false
    abs(integral_over_reference_element(f7, 3, order=7) - 1/720) >  1e-15 && return false
    
    failed = false
    try
        abs(integral_over_reference_element(f8, 3, order=8) - 1/990) >  1e-15
    catch e
        if isa(e, ErrorException) && occursin("Order 8", e.msg) && occursin("3D", e.msg)
            failed = true
        end
    end
    !failed && return false

    return true
end

@test test_quadrature()