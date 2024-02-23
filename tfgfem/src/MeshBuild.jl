
"""
    (nodes, elements, v_nodes, mesh) = MeshBuild(v_size, v_gap)

Given the size of the 3 dimensions (x,y,z) and the gap between points compute the
mesh.

The number of nodes and elements are also computed.

"""
function MeshBuild(v_size::Array{Int64, 1}, v_gap::Array{Float64,1})
    #ESste vector va a almacenar las coordenadas de los nodos
    v_nodes=[]
    #println("iniciando creacion de malla")
    """ DEFINIMOS EL TAMAÑO Y DIVISIONES DE NUESTRA ESTRUCTURA """
    #TAMAÑO
    sizex=v_size[1]
    sizey=v_size[2]
    sizez=v_size[3]

    #Tamaño de las divisiones. Vamos a imponer el tamaño deseado y a partir de él se
    # calcularan el numero de divisiones necesarias.
    gapx=v_gap[1]
    gapy=v_gap[2]
    gapz=v_gap[3]

    #DIVISIONES
    #=ndivx=5
    ndivy=5
    ndivz=20=#
    ndivx = round(Int64, sizex/gapx)
    ndivy = round(Int64, sizey/gapy)
    ndivz = round(Int64, sizez/gapz)
    #El numero de puntos es igual a (ndivx+1)*(ndivy+1)*(ndivz+1)

    #Esto es en caso de que queramos que el primer nodo esté desplazado del origen
    jumpx=0

    #=gapx=sizex/ndivx
    gapy=sizey/ndivy
    gapz=sizez/ndivz=#

    nodes = (ndivx+1)*(ndivy+1)*(ndivz+1)
    elements = ndivx*ndivy*ndivz*5 #Este código malla en tetraedros por lo que de cada cubo, sacamos 5 elementos.

    v_nodes = []
    for i in range(1,stop=ndivx+1)
        for j in range(1,stop=ndivy+1)
            for k in range(1,stop=ndivz+1)
                push!(v_nodes, [i*gapx+jumpx-sizex/2,-j*gapy,k*gapz])
            end
        end
    end


    mesh=[]
    #Este bucle recorre todos los elementos en direccion x. Los siguientes hacen lo mismo
    # en y y z.
    for i in range(1,stop=ndivx)

        for j in range(1,stop=ndivy)

            for k in range(1,stop=ndivz)
                global n_glob
                # numbering: i+j*nx+k*nx*ny
                # [i,j,k] es la posicion que ocupa el nodo local n(1,2,3...) en la malla
                n1=[i,j,k]
                n2=[i+1,j,k]
                n4=[i,j+1,k]
                n3=[i+1,j+1,k]
                n5=[i,j,k+1]
                n6=[i+1,j,k+1]
                n8=[i,j+1,k+1]
                n7=[i+1,j+1,k+1]

                n_loc=[n1,n2,n3,n4,n5,n6,n7,n8]
                n_glob=[]
                for l in range(1, stop=length(n_loc))
                    #n_glob.append(n_loc[l][0]*(ndivx+1)*(ndivy+1)+n_loc[l][1]*(ndivz+1)+n_loc[l][2])
                    push!(n_glob, (n_loc[l][1]-1)*(ndivz+1)*(ndivy+1)+(n_loc[l][2]-1)*(ndivz+1)+(n_loc[l][3]-1) + 1)
                end

                # El vec_aux nos va a permitir obtener los tetraedros a partir de los elementos cubicos.
                vec_aux=[[1,2,4,5],[2,3,4,7],[5,8,7,4],[5,7,6,2],[5,7,2,4]]
                for l in range(1,stop=length(vec_aux))
                    push!(mesh, [n_glob[vec_aux[l][1]], n_glob[vec_aux[l][2]], n_glob[vec_aux[l][3]], n_glob[vec_aux[l][4]]])
                end
            end
        end
    end

    npar=length(v_nodes)

    v_tetra=[]
    for i in range(1,stop=Int((npar)))

        for k in range(1,stop=3)
            push!(v_tetra, v_nodes[i][k])
        end
    end

    #fini= open(savedir*"\\"*"nodes.txt","w")
    # Check if mesh-data dir exists, create it otherwise
    if !isdir("mesh-data")
        mkdir("mesh-data")
    end
    fini = open("./mesh-data/nodes.txt", "w")
    for i in range(1,stop=length(v_nodes))

        for j in range(1,stop=length(v_nodes[i]))
            write(fini, string((v_nodes[i][j])))
            write(fini, "  ")
        end
        if i < length(v_nodes)
            write(fini, "\n")
        end
    end
    close(fini)



    # fini= open(savedir*"\\"*"mesh.txt","w")
    fini= open("./mesh-data/mesh.txt","w")

    mechanical_mesh = []
    for i in range(1,stop=length(mesh))
        for j in range(1,stop=length(mesh[i]))
            write(fini, string((mesh[i][j])))
            write(fini, "  ")
        end
        if i < length(mesh)
            write(fini, "\n")
        end

    end

    close(fini)

    return (nodes, elements, v_nodes, mesh)
end



"""
    (nodes, elements, v_nodes, mesh, mechanical_mesh) = MeshBuild(v_size, v_gap, E, ν)

Compute an array including mechanical properties such as Ypoung modulus and ν ratio

"""

function MeshBuild(v_size::Array{Float64, 1}, v_gap::Array{Float64,1}, E::Float64, ν::Float64)
    #ESste vector va a almacenar las coordenadas de los nodos
    v_nodes=[]
    #println("iniciando creacion de malla")
    """ DEFINIMOS EL TAMAÑO Y DIVISIONES DE NUESTRA ESTRUCTURA """
    #TAMAÑO
    sizex=v_size[1]
    sizey=v_size[2]
    sizez=v_size[3]

    #Tamaño de las divisiones. Vamos a imponer el tamaño deseado y a partir de él se
    # calcularan el numero de divisiones necesarias.
    gapx=v_gap[1]
    gapy=v_gap[2]
    gapz=v_gap[3]

    #DIVISIONES
    #=ndivx=5
    ndivy=5
    ndivz=20=#
    ndivx = round(Int64, sizex/gapx)
    ndivy = round(Int64, sizey/gapy)
    ndivz = round(Int64, sizez/gapz)
    #El numero de puntos es igual a (ndivx+1)*(ndivy+1)*(ndivz+1)

    #Esto es en caso de que queramos que el primer nodo esté desplazado del origen
    jumpx=0

    #=gapx=sizex/ndivx
    gapy=sizey/ndivy
    gapz=sizez/ndivz=#

    nodes = (ndivx+1)*(ndivy+1)*(ndivz+1)
    elements = ndivx*ndivy*ndivz*5 #Este código malla en tetraedros por lo que de cada cubo, sacamos 5 elementos.

    v_nodes = []
    for i in range(1,stop=ndivx+1)
        for j in range(1,stop=ndivy+1)
            for k in range(1,stop=ndivz+1)
                push!(v_nodes, [i*gapx+jumpx-sizex/2,-j*gapy,k*gapz])
            end
        end
    end


    mesh=[]
    #Este bucle recorre todos los elementos en direccion x. Los siguientes hacen lo mismo
    # en y y z.
    for i in range(1,stop=ndivx)

        for j in range(1,stop=ndivy)

            for k in range(1,stop=ndivz)
                global n_glob
                # numbering: i+j*nx+k*nx*ny
                # [i,j,k] es la posicion que ocupa el nodo local n(1,2,3...) en la malla
                n1=[i,j,k]
                n2=[i+1,j,k]
                n4=[i,j+1,k]
                n3=[i+1,j+1,k]
                n5=[i,j,k+1]
                n6=[i+1,j,k+1]
                n8=[i,j+1,k+1]
                n7=[i+1,j+1,k+1]

                n_loc=[n1,n2,n3,n4,n5,n6,n7,n8]
                n_glob=[]
                for l in range(1, stop=length(n_loc))
                    #n_glob.append(n_loc[l][0]*(ndivx+1)*(ndivy+1)+n_loc[l][1]*(ndivz+1)+n_loc[l][2])
                    push!(n_glob, (n_loc[l][1]-1)*(ndivz+1)*(ndivy+1)+(n_loc[l][2]-1)*(ndivz+1)+(n_loc[l][3]-1) + 1)
                end

                # El vec_aux nos va a permitir obtener los tetraedros a partir de los elementos cubicos.
                vec_aux=[[1,2,4,5],[2,3,4,7],[5,8,7,4],[5,7,6,2],[5,7,2,4]]
                for l in range(1,stop=length(vec_aux))
                    push!(mesh, [n_glob[vec_aux[l][1]], n_glob[vec_aux[l][2]], n_glob[vec_aux[l][3]], n_glob[vec_aux[l][4]]])
                end
            end
        end
    end

    npar=length(v_nodes)

    v_tetra=[]
    for i in range(1,stop=Int((npar)))

        for k in range(1,stop=3)
            push!(v_tetra, v_nodes[i][k])
        end
    end

    # Check if mesh-data dir exists, create it otherwise
    if !isdir("mesh-data")
        mkdir("mesh-data")
    end
    #fini= open(savedir*"\\"*"nodes.txt","w")
    fini = open("./mesh-data/nodes.txt", "w")
    for i in range(1,stop=length(v_nodes))

        for j in range(1,stop=length(v_nodes[i]))
            write(fini, string((v_nodes[i][j])))
            write(fini, "  ")
        end
        if i < length(v_nodes)
            write(fini, "\n")
        end
    end
    close(fini)


    #ecrit=0.000134
    #efin=0.00015
    #subdivided=0

    # fini= open(savedir*"\\"*"mesh.txt","w")
    fini= open("./mesh-data/mesh.txt","w")
    # fmech= open(savedir*"\\"*"mechanical_mat.txt","w"
    fmech= open("./mesh-data/mechanical_mat.txt","w")
    mechanical_mesh = []

    for i in range(1,stop=length(mesh))
        for j in range(1,stop=length(mesh[i]))
            write(fini, string((mesh[i][j])))
            write(fini, "  ")
        end
        if i < length(mesh)
            write(fini, "\n")
        end

        #push!(mechanical_mesh, [E, ν, ecrit, efin, subdivided])
        push!(mechanical_mesh, [E, ν])
        write(fmech, string( E  ))
        write(fmech, "  ")
        write(fmech, string( ν  ))
        #=write(fmech, "  ")
        write(fmech, string( ecrit  ))
        write(fmech, "  ")
        write(fmech, string( efin  ))
        write(fmech, "  ")
        write(fmech, string( subdivided  ))
        =#
        write(fmech, "\n")
    end

    close(fini)
    close(fmech)

    return (nodes, elements, v_nodes, mesh, mechanical_mesh)
end


"""
    (E_mesh, ν_mesh) = MechanicalMesh(v_mechanical)

Separate the different mechanical properties included in the argument array
"""
#Funcion para inicializar las propiedades mecanicas de los elementos de la malla

function MechanicalMesh(v_mechanical)
    elements = length(v_mechanical)

    E_mesh=[]
    for i in range(1,stop=elements)
        E1=v_mechanical[i][1]
        push!(E_mesh, E1)
    end

    ν_mesh=[]
    for i in range(1,stop=elements)
        ν1=v_mechanical[i][2]
        push!(ν_mesh, ν1)
    end

    return (E_mesh, ν_mesh)
end


"""
    E_nodes = NodalYoungModule(nodal_connectivity, E_mesh)

Compute the equivalent Young modulus of each node using the modulus of the elements
surrounding it.
It can be used for further printing in programs such as Paraview.

"""
function NodalYoungModule(nodal_connectivity, E_mesh)

    nodes = length(nodal_connectivity)
    nodal_E = []

    # recorremos cada nodo
    for i in range(1, stop=nodes)
        Ei = 0
        for j in range(1, stop=length(nodal_connectivity[i]))
            Ei += E_mesh[nodal_connectivity[i][j]] / (length(nodal_connectivity[i]))

        end
        push!(nodal_E, Ei)
    end

    return nodal_E

end
