
function CubeLayeringZ(v_nodes::Array{Any,1}, v_mesh::Array{Any,1},
    mesh_distribution::Array{Any,1})

    num_layers = length(mesh_distribution[3]) - 1
    layers = []
    v_mesh_layer = []
    for ilayer = 1:(num_layers)

        z1 = mesh_distribution[3][ilayer + 1]
        z2 = mesh_distribution[3][ilayer]

        nodes_z1 = NodeLayerZ(z1, v_nodes)
        nodes_z2 = NodeLayerZ(z2, v_nodes)

        (element_layer, mesh_layer) = ElementLayer(v_nodes, v_mesh, nodes_z1, nodes_z2)

        push!(layers, element_layer)
        push!(v_mesh_layer, mesh_layer)

    end



    return layers, v_mesh_layer #v_layer_nodes

end


##
function ElementLayer(v_nodes::Array{Any,1}, v_mesh::Array{Any,1},
    node_layer1::Array{Any,1}, node_layer2::Array{Any,1})

    elements = length(v_mesh)

    element_layer = []
    v_mesh_layer = []
    for i = 1:elements

        check1 = 0
        check2 = 0
        for n = 1:4

            node =v_mesh[i][n]
            if count(i->(i==node), node_layer1) != 0
                check1 = 1
            end

            if count(i->(i==node), node_layer2) != 0
                check2 = 1
            end

        end

        if check1==1 && check2==1
            push!(element_layer, i)
            push!(v_mesh_layer, v_mesh[i])
        end

    end


    return element_layer, v_mesh_layer

end
