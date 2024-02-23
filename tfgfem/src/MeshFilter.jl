function MeshFilter(nodal_connectivity, elements)
    """Since nodal_connectivity takes into account the elements per node, this function filters out
    some given *elements* from the corresponding nodes, such that the final Young Modulus
    smoothing is performed only in working elements
    """
    for node in nodal_connectivity
        filter!(e -> e ∈ elements, node)
    end
    return nodal_connectivity
end
