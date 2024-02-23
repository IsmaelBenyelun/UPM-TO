
function RestrictedDOF(BCMatrix::Array{Any, 1}, u::Array{Float64,1})

    restr_nodes = length(BCMatrix)
    restricted_DOF = []

    for i=1:restr_nodes

        restr_node = BCMatrix[i][1]
        x_displ = BCMatrix[i][5]
        y_displ = BCMatrix[i][6]
        z_displ = BCMatrix[i][7]


        if BCMatrix[i][2] == 1
            u[Int(restr_node*3-2)] = x_displ
            push!(restricted_DOF, Int(restr_node*3-2))
        end

        if BCMatrix[i][3] == 1
            u[Int(restr_node*3-1)] = y_displ
            push!(restricted_DOF, Int(restr_node*3-1))
        end

        if BCMatrix[i][4] == 1
            u[Int(restr_node*3)] = z_displ
            push!(restricted_DOF, Int(restr_node*3))
        end

    end


    return restricted_DOF, u

end
