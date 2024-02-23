"""
    (Kglob, B) = Kglobal(dof_node, v_nodes, v_mesh, E_mesh, pois_mesh)

Compute the global K matrix of a structure and the matrix B
dof_node: degrees of freedom per node
E_mesh and pois_mesh: arrays including the mesh mechanical properties
"""

using SparseArrays, Base.Threads

function Kglobal(dof_node::Int64, v_nodes::Array{Any,1}, v_mesh::Array{Any,1},
    E_mesh::Array{Any,1}, pois_mesh::Array{Any,1})

    nodes = length(v_nodes)
    elements = length(v_mesh)
    dof_node = 3
    ndof = nodes*dof_node
    Kglob = spzeros(ndof,ndof)
    collect_aux_vec_coord = []

    for l in range(1, stop=elements)
        global B
        el_loc=l
        E_el=E_mesh[l]    #Young module for element l
        pois_el =  pois_mesh[l]  #poisson for element l
        aux_vec_coord=[]

        # VECTOR CON LAS COORDENADAS LOCALES DEL ELEMENTO TETRAHEDRICO
        #aux_vec_coord=[[]]
        #=Este for recorre en las cuatro columnas que tiene el array v_mesh para
        ir almcenando los nodos correspondientes al elemento el_loc con el que estamos
        en cada momento.
        En el segundo for se recogen las coordenadas de dicho nodo y se almacenan
        en el vector auxiliar aux_vec_coord=#
        for j2 in range(1,stop=4)
            #En naux se almacena
            naux=v_mesh[el_loc][j2] #NODO DE LA CONECTIVIDAD
            for k2 in range(1,stop=3)  #recorre las 3 columnas de v_nodes
                push!(aux_vec_coord, v_nodes[naux][k2])
            end
        end
        #println("auxilary coodinate vector", aux_vec_coord)
        push!(collect_aux_vec_coord, aux_vec_coord)
        # FUNCION EXTERdNA PARA CALCULAR LA MATRIZ DE RIGIDEZ LOCAL DEL ELEMENTO
        (Kloc, B, voltet, D, straMat) = KlocTetra(aux_vec_coord, E_el, pois_el)

        #=Vamos a definir ahora un vector que recoja los grados de libertad correspondientes
        a cada elemento el_loc del bucle en el que estamos. Va a tener 12 terminos,
        ya que cada elemento tiene 12 gdl (4elementos x 3 gdl/nodo)=#
        vec_loc_ind=[]
        #=j2 selecciona uno de los 4 nodos del elemento el_loc. Esto lo hace buscando
        en la fila correspondiente de la matriz v_mesh=#
        for j2 in range(1,stop=4)
            #j3 selecciona el grado de libertad del nodo deseado
            for j3 = (2:-1:0)
                push!(vec_loc_ind, v_mesh[el_loc][j2]*3-j3)
            end
        end

        """ENSAMBLAJE DE LA LOCAL EN LA GLOBAL"""
        #=los contadores iloc y jloc van a recorrer la Kloc a la vez que van situando
        los valores de la Kloc en la correspondiente posicion en la Kglob. Para ello,
        utiliza el vector vec_loc_ind calculado anteriormente=#
        for iloc in range(1,stop=length(Kloc[1,:])) # ensamblaje
        #for iloc in 1:size(Kloc)[1]
            for jloc in range(1,stop=length(Kloc[:, iloc]))
            #for jloc in 1:size(Kloc)[1]

                Kglob[vec_loc_ind[iloc], vec_loc_ind[jloc]] += voltet*Kloc[iloc, jloc]

            end
        end

    end

    return (Kglob, B, collect_aux_vec_coord)
end


"""
    (KglobSp, B) = KglobalSp(dof_node, v_nodes, v_mesh, E_mesh, pois_mesh)

Compute the global K matrix of a structure and the matrix B
dof_node: degrees of freedom per node
E_mesh and pois_mesh: arrays including the mesh mechanical properties
"""

function KglobalSp(dof_node::Int64, v_nodes::Array{Any,1}, v_mesh::Array{Any,1},
   E_mesh::Array{Any,1}, pois_mesh::Array{Any,1})

    nodes = length(v_nodes)
    elements = length(v_mesh)
    dof_node = 3
    ndof = nodes*dof_node
    collect_aux_vec_coord = []

    I = []
    J = []
    V = []
    for l in range(1, stop=elements)
        global B
        el_loc=l
        E_el=E_mesh[l]    #Young module for element l
        pois_el =  pois_mesh[l]  #poisson for element l
        aux_vec_coord=[]

        # VECTOR CON LAS COORDENADAS LOCALES DEL ELEMENTO TETRAHEDRICO
        #aux_vec_coord=[[]]
        #=Este for recorre en las cuatro columnas que tiene el array v_mesh para
        ir almcenando los nodos correspondientes al elemento el_loc con el que estamos
        en cada momento.
        En el segundo for se recogen las coordenadas de dicho nodo y se almacenan
        en el vector auxiliar aux_vec_coord=#
        for j2 in range(1,stop=4)
            #En naux se almacena
            naux=v_mesh[el_loc][j2]
            for k2 in range(1,stop=3)  #recorre las 3 columnas de v_nodes
                push!(aux_vec_coord, v_nodes[naux][k2])
            end
        end
        #println("auxilary coodinate vector", aux_vec_coord)
        push!(collect_aux_vec_coord, aux_vec_coord)
        # FUNCION EXTERdNA PARA CALCULAR LA MATRIZ DE RIGIDEZ LOCAL DEL ELEMENTO
        (Kloc, B, voltet, D, straMat) = KlocTetra(aux_vec_coord, E_el, pois_el)

        #=Vamos a definir ahora un vector que recoja los grados de libertad correspondientes
        a cada elemento el_loc del bucle en el que estamos. Va a tener 12 terminos,
        ya que cada elemento tiene 12 gdl (4elementos x 3 gdl/nodo)=#
        vec_loc_ind=[]
        #=j2 selecciona uno de los 4 nodos del elemento el_loc. Esto lo hace buscando
        en la fila correspondiente de la matriz v_mesh=#
        for j2 in range(1,stop=4)
            #j3 selecciona el grado de libertad del nodo deseado
            for j3 = (2:-1:0)
                push!(vec_loc_ind, v_mesh[el_loc][j2]*3-j3)
            end
        end

        """ENSAMBLAJE DE LA LOCAL EN LA GLOBAL"""
        #=los contadores iloc y jloc van a recorrer la Kloc a la vez que van situando
        los valores de la Kloc en la correspondiente posicion en la Kglob. Para ello,
        utiliza el vector vec_loc_ind calculado anteriormente=#

        for iloc in range(1,stop=length(Kloc[1,:])) # ensamblaje
        #for iloc in 1:size(Kloc)[1]
            for jloc in range(1,stop=length(Kloc[:, iloc]))
            #for jloc in 1:size(Kloc)[1]
                value = voltet*Kloc[iloc, jloc]

                if value != 0
                    push!(I, vec_loc_ind[iloc])
                    push!(J, vec_loc_ind[jloc])
                    push!(V, value)
                end

            end
        end

    end

    KglobalSparse = sparse(I,J,V,ndof,ndof)

    return (KglobalSparse, B)
end
