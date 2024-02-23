

function KglobalLL(Kglob, u, restrict_DOF)

    # u array contains de pre-included displacements
    u = convert(Array{Float64,1}, u)
    u_mat=copy(u)

    #num_FREE = number of free degrees of  freedom
    #num_RESTR = number of restricted degrees of freedom
    free_DOF = setdiff(1:132, restrict_DOF)

    num_FREE = length(free_DOF)
    num_RESTR = length(restrict_DOF)

    Kglob_FF = zeros(num_FREE, num_FREE)
    Kglob_FR = zeros(num_FREE, num_RESTR)
    Kglob_RR = zeros(num_RESTR, num_RESTR)

    # Force vector including free degrees of freedom.
    # If nothing is said, they are zero.
    f_F = zeros(num_FREE)

    # Force vector including restricted degrees of freedom.
    # They are unknown.
    f_R = zeros(num_RESTR)

    # Building the Kglob_FF matrix.
    for irow = 1:num_FREE
        row = free_DOF[irow]
        for icolumn = 1:num_FREE
            column = free_DOF[icolumn]
            Kglob_FF[irow, icolumn] = Kglob[row, column]
        end
    end

    # Building the Kglob_RR matrix.
    for irow = 1:num_RESTR
        row = restrict_DOF[irow]
        for icolumn = 1:num_RESTR
            column = restrict_DOF[icolumn]
            Kglob_RR[irow, icolumn] = Kglob[row, column]
        end
    end

    # Building the Kglob_FR matrix.
    for irow = 1:num_FREE
        row = free_DOF[irow]
        for icolumn = 1:num_RESTR
            column = restrict_DOF[icolumn]
            Kglob_FR[irow, icolumn] = Kglob[row, column]
        end
    end
    #Kglob_RF = transpose(Kglob_FR)


    u_R = zeros(num_RESTR)
    for irow = 1:num_RESTR
        row = restrict_DOF[irow]
        u_R[irow] = u[row]
    end

    return (Kglob_FF, Kglob_FR, Kglob_RR, u_R)

end
