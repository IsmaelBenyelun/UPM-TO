"""
    circle = CircleZ(v_nodes, center_coords, radius)

Select a group of nodes inside a circle at a given Z-height

center_coords: mesh coordinates of the center of the circle
"""

function CircleZ(v_nodes::Array{Any,1}, center_coords::Array{Float64,1}, radius::Float64)

    z_layer = center_coords[3]
    center_node = NodeSelector(v_nodes, center_coords)

    # Computing the layer where the circle is going to be selected
    nodes = length(v_nodes)
    layer_nodes = []

    for i=1:nodes

        z_i = v_nodes[i][3]

        if z_i == z_layer
            push!(layer_nodes, i)
        end

    end

    # Choosing the nodes inside the selected radius
    x_c = v_nodes[center_node][1]
    y_c = v_nodes[center_node][2]

    circle = []

    for i=1:length(layer_nodes)
        i_node = layer_nodes[i]
        x_i = v_nodes[i_node][1]
        y_i = v_nodes[i_node][2]
        # Computing distance between reference node and node i
        distance = sqrt((x_c - x_i)^2 + (y_c - y_i)^2)

        if distance <= radius
            push!(circle, i_node)
        end

    end

    return circle

end

"""
    circle = Circle(v_nodes, indicator, center_coords, radius)

Select a group of nodes inside a circle at a given layer. This layer can be x-constant,
y-constant or z-constant. This information is provided by the indicator vector

center_coords: mesh coordinates of the center of the circle
"""

function Circle(v_nodes::Array{Any,1}, indicator::Array{Int64,1},
    center_coords::Array{Float64,1}, radius::Float64)

    if indicator[1] == 1 #x-constant
        layer = center_coords[1]
    elseif indicator[2] == 1 #y-constant
        layer = center_coords[2]
    elseif indicator[3] == 1 #z-constant
        layer = center_coords[3]
    end

    center_node = NodeSelector(v_nodes, center_coords)

    # Computing the layer where the circle is going to be selected
    nodes = length(v_nodes)
    layer_nodes = []

    for i=1:nodes

        x_i = v_nodes[i][1]
        y_i = v_nodes[i][2]
        z_i = v_nodes[i][3]

        if (indicator[1] == 1) && (x_i == layer)
            push!(layer_nodes, i)
        elseif (indicator[2] == 1) && (y_i == layer)
            push!(layer_nodes, i)
        elseif (indicator[3] == 1) && (z_i == layer)
            push!(layer_nodes, i)
        end

    end

    # Choosing the nodes inside the selected radius
    x_c = v_nodes[center_node][1]
    y_c = v_nodes[center_node][2]
    z_c = v_nodes[center_node][3]

    circle = []

    for i=1:length(layer_nodes)
        i_node = layer_nodes[i]

        if (indicator[1] == 1)
            z_i = v_nodes[i_node][3]
            y_i = v_nodes[i_node][2]
            # Computing distance between reference node and node i
            distance = sqrt((z_c - z_i)^2 + (y_c - y_i)^2)

        elseif (indicator[2] == 1)
            x_i = v_nodes[i_node][1]
            z_i = v_nodes[i_node][3]
            # Computing distance between reference node and node i
            distance = sqrt((x_c - x_i)^2 + (z_c - z_i)^2)

        elseif (indicator[3] == 1)
            x_i = v_nodes[i_node][1]
            y_i = v_nodes[i_node][2]
            # Computing distance between reference node and node i
            distance = sqrt((x_c - x_i)^2 + (y_c - y_i)^2)

        end

        if distance <= radius
            push!(circle, i_node)
        end

    end

    return circle

end


"""
    circunference = CircunferenceZ(v_nodes, center_coords, rmin, rmax)

"""
function CircunferenceZ(v_nodes::Array{Any,1}, center_coords::Array{Float64,1},
    rmin::Float64, rmax::Float64)

    z_layer = center_coords[3]
    center_node = NodeSelector(v_nodes, center_coords)

    # Computing the layer where the circle is going to be selected
    nodes = length(v_nodes)
    layer_nodes = []

    for i=1:nodes

        z_i = v_nodes[i][3]
        if z_i == z_layer
            push!(layer_nodes, i)
        end

    end

    # Choosing the nodes between the two radius
    x_c = v_nodes[center_node][1]
    y_c = v_nodes[center_node][2]

    circunference = []

    for i=1:length(layer_nodes)
        i_node = layer_nodes[i]
        x_i = v_nodes[i_node][1]
        y_i = v_nodes[i_node][2]
        # Computing distance between reference node and node i
        distance = sqrt((x_c - x_i)^2 + (y_c - y_i)^2)

        if (rmin < distance < rmax)
            push!(circunference, i_node)
        end

    end

    return circunference

end


"""
    select_node = NodeSelector(v_nodes, node_coords)

Select a node from the mesh with the coordinates specified.

"""
function NodeSelector(v_nodes, node_coords)

    nodes = length(v_nodes)
    select_node = 0

        # Setting the selected coordinates
        x_s = node_coords[1]
        y_s = node_coords[2]
        z_s = node_coords[3]

        for i = 1:nodes

            x_i = v_nodes[i][1]
            y_i = v_nodes[i][2]
            z_i = v_nodes[i][3]

            if (x_i == x_s) && (y_i == y_s) && (z_i == z_s)
                select_node = i
            end

        end


        if select_node == 0
            println("no node was found on those coordinates")
        end

    return select_node

end

"""
    node_group = NodeGroup(v_nodes, dimension, coords)

Select a group of nodes that meet certain coordinates.
dimension: 3-component vector to activate(1) or desactivate(2) dimensions
x, y or z. If desactivated, that coordinate won't be taken into account.

"""

function NodeGroup(v_nodes::Vector{Any}, dimension::Array{Int64,1},
    coords::Array{Float64})

    nodes = length(v_nodes)

    dimx = dimension[1]
    dimy = dimension[2]
    dimz = dimension[3]

    node_group = []
    for i = 1:nodes

        xloc=v_nodes[i][1]
        yloc=v_nodes[i][2]
        zloc=v_nodes[i][3]

        if (dimx==1 && dimy==0 && dimz==0) && (xloc==coords[1])
            push!(node_group, i)

        elseif (dimx==0 && dimy==1 && dimz==0) && (yloc==coords[2])
            push!(node_group, i)

        elseif (dimx==0 && dimy==0 && dimz==1) && (zloc==coords[3])
            push!(node_group, i)

        elseif (dimx==1 && dimy==1 && dimz==0) && (xloc==coords[1] && yloc==coords[2])
            push!(node_group, i)

        elseif (dimx==0 && dimy==1 && dimz==1) && (yloc==coords[2] && zloc==coords[3])
            push!(node_group, i)

        elseif (dimx==1 && dimy==0 && dimz==1) && (xloc==coords[1] && zloc==coords[3])
            push!(node_group, i)

        elseif (dimx==1 && dimy==1 && dimz==1) && (xloc==coords[1] && yloc == coords[2] && zloc==coords[3])
            push!(node_group, i)

        end

    end


    return node_group

end

"""
    node_layer = NodeLayerZ(z_coord, v_nodes)

Select a layer of nodes located at the same z-coordinate
"""

function NodeLayerZ(z_coord::Float64, v_nodes)

    nodes = length(v_nodes)
    node_layer = []
    for i = 1:nodes

        z_node = v_nodes[i][3]
        if z_node == z_coord
            push!(node_layer, i)
        end

    end

    return node_layer

end
