
function invDrMatrix(ν::Float64)

    d1 = 1
    d2 = -ν
    d3 = 2*(1 + ν)
    # d3 = (1 + ν) # Wrong!

    # Compute D
    invDr = [
        d1 d2 d2 0 0 0;
        d2 d1 d2 0 0 0;
        d2 d2 d1 0 0 0;
        0 0 0 d3 0 0;
        0 0 0 0 d3 0;
        0 0 0 0 0 d3
    ]
    return invDr
end
