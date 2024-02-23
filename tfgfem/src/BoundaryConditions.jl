#=Funcion para aplicar las condiciones de contorno.
Se verifican qué nodos cumplen las condiciones impuestas y se almacenan en los vectores
de salida de la funcion.
Pueden ser condiciones de contorno de fuerzas o esplazamientos. Nosotros las
vamos a aplicar de desplazamiento
=#

function BoundaryConditions(v_nodes::Vector{Any}, elements::Int64,
    pos_BC1::Array{Float64,1}, pos_BC2::Array{Float64,1}, BC1::Float64, BC2::Float64)

    nodes = length(v_nodes)

    for l in range(1, stop=elements)
        global nodes_BC1, nodes_BC2, u, restrict_DOF
        #Nodos primera BC
        nodes_BC1=[]
        #Nodos inferiores
        nodes_BC2=[]
        #Displacements vector
        u=[]
        #=Vector encargado de eliminar las filas y columnas asociadas a las
        BC en la KGlob=#
        restrict_DOF=[]


        for i in range(1,stop=nodes)
            #println("nodo", i)
            xloc=v_nodes[i][1]
            yloc=v_nodes[i][2]
            zloc=v_nodes[i][3]

            if xloc==pos_BC2[1] && zloc==pos_BC2[3]
                #=Si el nodo cumple las condiciones de posicion anteriores, lo guardo en el
                vector de nodos inferiores.=#
                #println("hola", restrict_DOF)
                push!(nodes_BC2, i)
                #println("nodosBC2", nodes_BC2)
                #= k counter will add the restricted degrees of freeedom from
                node i to the vector restrict_DOF
                Depending which DOF we want to restrict, we must change the number
                substracting
                -2 : x displacement
                -1 : y displacement
                0 : z displacement
                =#
                #for k=(2:-1:0)
                    push!(restrict_DOF, i*3)
                #end
            end
            #println("remove", restrict_DOF)


            if xloc==pos_BC1[1]
                #=Si el nodo cumple las condiciones de posicion anteriores, lo guardo en el
                vector de nodos superiores.=#
                push!(nodes_BC1, i)
                #println("nodosBC1", nodes_BC1)
                # Embedded condition: the 3 DOF are restricted (for loop)
                for k=(2:-1:0)
                    #=La funcion count cuenta el numero de veces que ij==i*3+k. Es decir,
                    este contador se asegura que este gdl global no se haya añadido ya
                    al vectore restrict_DOF=#
                    if count(ij->(ij==i*3-k), restrict_DOF)==0
                        push!(restrict_DOF, i*3-k)
                    end
                end
            end

            #println("hola", restrict_DOF)

            #Initializing both vectors
            for j in range(1,stop=3)
                push!(u, 0.0)
                #push!(dx, 0.0)
            end
        end

        #=Application of BC1=#
        for i in range(1,stop=length(nodes_BC1))

            iloc1=nodes_BC1[i]
            #= In the correspondent position (node iloc1), y displacement BC is
            applied=#
            u[iloc1*3-1]= BC1
        end

        # Application of BC2
        for i in range(1,stop=length(nodes_BC2))

            iloc2=nodes_BC2[i]
            #=En la posicion correspondiente (número del nodo, iloc2), se aplica
            una condicion de contorno de desplazamiento vertical=#
            u[iloc2*3]= BC2
        end
    end

    return (nodes_BC1, nodes_BC2, restrict_DOF, u)

end


"""
    BCMat = BCMatrix(BCMat, nodes_BC, condition, BC)

Compute a matrix for applying the boundary conditions. This matrix must be an input
(BCMat)

nodes_BC: nodes where the BC will be applied
condition: 1 and 0 arrays that select which displacements will be affected by the BC
BC: the value of the selected displacements

"""
function BCMatrix(BCMat::Array{Any, 1}, nodes_BC::Array{Any,1},
    condition::Array{Int64,1}, BC::Array{Float64,1})

    for i=1:length(nodes_BC)

        restr_node = nodes_BC[i]

        # x displacement condition
        if condition[1]==1
            restr_x = 1
            BCx = BC[1]
        else
            restr_x = 0
            BCx = 0.0
        end

        # y displacement
        if condition[2]==1
            restr_y = 1
            BCy = BC[2]
        else
            restr_y = 0
            BCy = 0.0
        end

        # z displacement
        if condition[3]==1
            restr_z = 1
            BCz = BC[3]
        else
            restr_z = 0
            BCz = 0.0
        end

        restr_node = Int(restr_node)
        push!(BCMat, [restr_node, restr_x, restr_y, restr_z, BCx, BCy, BCz])
    end

    return BCMat

end

function BCMatrix(BCMat::Array{Any, 1}, nodes_BC::Array{Any, 1},
    condition::Array{Int64, 1}, BC::Array{Array{Float64, 1}, 1})

    for i=1:length(nodes_BC)

        restr_node = nodes_BC[i]

        # x displacement condition
        if condition[1]==1
            restr_x = 1
            BCx = BC[i][1]
        else
            restr_x = 0
            BCx = 0.0
        end

        # y displacement
        if condition[2]==1
            restr_y = 1
            BCy = BC[i][2]
        else
            restr_y = 0
            BCy = 0.0
        end

        # z displacement
        if condition[3]==1
            restr_z = 1
            BCz = BC[i][3]
        else
            restr_z = 0
            BCz = 0.0
        end

        restr_node = Int(restr_node)
        push!(BCMat, [restr_node, restr_x, restr_y, restr_z, BCx, BCy, BCz])
    end

    return BCMat

end
