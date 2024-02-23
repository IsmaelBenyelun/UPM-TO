"""UNA VEZ TENEMOS LA MATRIZ DE RIGIDEZ Y LAS CONDICIONES DE CONTORNO RESOLVEMOS"""

function SolveFEM(Kglob, u, f, restrict_DOF)
    num_DOF = length(Kglob[:, 1])
    free_DOF = setdiff(1:num_DOF, restrict_DOF)

    # Force vector including free degrees of freedom.
    f_F = f[free_DOF]

    # Assign the Kglob_FF, Kglob_RR and Kglob_FR matricesm and u_R vector
    Kglob_FF = Kglob[free_DOF, free_DOF]
    Kglob_RR = Kglob[restrict_DOF, restrict_DOF]
    Kglob_FR = Kglob[free_DOF, restrict_DOF]

    u_R = u[restrict_DOF]

    # Solve
    b = f_F - Kglob_FR*u_R
    # u_F = Kglob_FF \ b #it is basically the same
    u_F = cholesky(Kglob_FF) \ b
    # Assign
    u[free_DOF] = u_F
    f_R = transpose(Kglob_FR) * u_F + Kglob_RR * u_R
    return (Kglob_FF, Kglob_FR, Kglob_RR, u_F, u_R, u, f_F, f_R)

end
