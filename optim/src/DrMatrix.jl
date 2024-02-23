
function DrMatrix(ν::Float64)

    d1 = (1 - ν) / ((1 - 2*ν) * (1 + ν))
    d2 = ν / ((1 - 2*ν)*(1 + ν))
    d3 = 1 / (2*(1 + ν))
    # d3 = 1 / (1 + ν) # Wrong! Voigt notation involves γ_xy, γ_xz and γ_yz (not ε_ij)

    # Compute D
    Dr = [
        d1 d2 d2 0 0 0;
        d2 d1 d2 0 0 0;
        d2 d2 d1 0 0 0;
        0 0 0 d3 0 0;
        0 0 0 0 d3 0;
        0 0 0 0 0 d3
    ]
    return Dr
end
