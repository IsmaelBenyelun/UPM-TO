"""
    (element_info, connect_element_glob) = ElementConnectivity(v_nodes, v_mesh)

Element_info::{Array{Array{Array}}}. The second array corresponds to an element.
Inside, there are four vectors than include the coordinates of each element's node

"""

function ElementConnectivity(v_nodes, v_mesh)

    element_info=[]
    connect_element_glob=[]

    elements = length(v_mesh)
    for i in range(1, stop=elements)

        push!(element_info, [])

        push!(connect_element_glob,v_mesh[i][1])
        push!(element_info[i],[v_nodes[v_mesh[i][1]][1], v_nodes[v_mesh[i][1]][2], v_nodes[v_mesh[i][1]][3]])
        #push!(element_info,v_nodes[v_mesh[i][1]][1])
        #push!(element_info,v_nodes[v_mesh[i][1]][2])
        #push!(element_info,v_nodes[v_mesh[i][1]][3])

        push!(connect_element_glob,v_mesh[i][2])
        push!(element_info[i],[v_nodes[v_mesh[i][2]][1], v_nodes[v_mesh[i][2]][2], v_nodes[v_mesh[i][2]][3]])
        #push!(element_info,v_nodes[v_mesh[i][2]][1])
        #push!(element_info,v_nodes[v_mesh[i][2]][2])
        #push!(element_info,v_nodes[v_mesh[i][2]][3])

        push!(connect_element_glob,v_mesh[i][3])
        push!(element_info[i],[v_nodes[v_mesh[i][3]][1], v_nodes[v_mesh[i][3]][2], v_nodes[v_mesh[i][3]][3]])
        #push!(element_info,v_nodes[v_mesh[i][3]][1])
        #push!(element_info,v_nodes[v_mesh[i][3]][2])
        #push!(element_info,v_nodes[v_mesh[i][3]][3])

        push!(connect_element_glob,v_mesh[i][4])
        push!(element_info[i],[v_nodes[v_mesh[i][4]][1], v_nodes[v_mesh[i][4]][2], v_nodes[v_mesh[i][4]][3]])
        #push!(element_info,v_nodes[v_mesh[i][4]][1])
        #push!(element_info,v_nodes[v_mesh[i][4]][2])
        #push!(element_info,v_nodes[v_mesh[i][4]][3])
    end
    return (element_info, connect_element_glob)
end


"""
    element_center_coord = ElementCenterCoord(element_info)

Compute the coordinates of the centes of the elements of a mesh
"""

function ElementCenterCoord(element_info)

    element_center_coord = []

    elements = length(element_info)

    for i=1:elements

        nodes_el = length(element_info[i])

        # X coordinate
        x_el = 0.0
        for j=1:nodes_el
            x_el = x_el + (element_info[i][j][1]) / 4.0
        end

        # Y coordinste
        y_el = 0.0
        for j=1:nodes_el
            y_el = y_el + (element_info[i][j][2]) / 4.0
        end

        # Z coordinate
        z_el = 0.0
        for j=1:nodes_el
            z_el = z_el + (element_info[i][j][3]) / 4.0
        end

        push!(element_center_coord, [x_el, y_el, z_el])

    end

    return element_center_coord

end



##
"""
    nodal_connectivity = NodalConnectivity(v_nodes, v_mesh)

Compute the mesh connectivity. For each node, the functions recognizes which
elements include that node

"""

function NodalConnectivity(v_nodes, v_mesh)

    nodal_connectivity=[]

    for i in range(1,stop=length(v_nodes))
        push!(nodal_connectivity, [])
    end
    # i recorre desde 1 hasta el numero de elementos
    for i in range(1,stop=length(v_mesh))
        # j va a recorrer las columnas de v_mesh
        for j in range(1,stop=length(v_mesh[1]))
            #println(nod_connectivity[v_mesh[i,j]])
            push!(nodal_connectivity[v_mesh[i][j]], i)
        end
    end

    return (nodal_connectivity)

end


"""
    coords_element = CoordsElement(v_mesh, v_nodes, el_loc)

Compute the coordinates of the element number el_loc of the mesh
"""

function CoordsElement(v_mesh::Array{Any,1}, v_nodes::Array{Any,1}, el_loc::Int64)

    coords_elem = []
    for j2 in range(1,stop=4)

        node=v_mesh[el_loc][j2]
        for k2 in range(1,stop=3)
            push!(coords_elem, v_nodes[node][k2])
        end
    end


    return coords_elem

end
