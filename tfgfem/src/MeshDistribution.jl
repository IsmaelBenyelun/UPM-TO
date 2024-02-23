"""
    mesh_distribution = MeshDistribution(v_size, v_gap, jumpx)

Compute where nodes are situated in each direction (x,y,z)
"""
function MeshDistribution(v_size, v_gap, jumpx)

    sizex=v_size[1]
    sizey=v_size[2]
    sizez=v_size[3]

    gapx=v_gap[1]
    gapy=v_gap[2]
    gapz=v_gap[3]

    #DIVISIONS
    ndivx = round(Int64, sizex/gapx)
    ndivy = round(Int64, sizey/gapy)
    ndivz = round(Int64, sizez/gapz)

    mesh_distribution = []

    push!(mesh_distribution, [])
    for i in range(1,stop=ndivx+1)
        push!(mesh_distribution[1], i*gapx+jumpx-sizex/2)
    end

    push!(mesh_distribution, [])
    for j in range(1,stop=ndivy+1)
        push!(mesh_distribution[2], -j*gapy)
    end

    push!(mesh_distribution, [])
    for k in range(1,stop=ndivz+1)
        push!(mesh_distribution[3], k*gapz)
    end


    return mesh_distribution

end
