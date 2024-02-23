
function FEAnalysis(v_nodes::Array{Any,1}, v_mesh::Array{Any,1}, E_mesh::Array{Any,1},
    pois_mesh::Array{Any,1}, u::Array{Float64,1}, f::Array{Float64,1},
    restricted_DOF::Array{Any,1}; print_f_R::Bool=false)

    # Obtaining the nodal connectivity
    nodal_connectivity = NodalConnectivity(v_nodes, v_mesh)

    # Obtaining the elemental connectivity
    element_info, connect_element_glob = ElementConnectivity(v_nodes, v_mesh)

    # Global K Matrix
    dof_node = 3
    # print("Assembly..."); @time Kglob, B = Kglobal(dof_node, v_nodes, v_mesh, E_mesh, pois_mesh)
    Kglob, B = Kglobal(dof_node, v_nodes, v_mesh, E_mesh, pois_mesh)
    num_DOF = length(Kglob[:, 1])

    # print("Solve FEM..."); @time (Kglob_FF, Kglob_FR, Kglob_RR, u_F, u_R, u, f_F, f_R) = SolveFEM(Kglob, u, f, restricted_DOF)
    (Kglob_FF, Kglob_FR, Kglob_RR, u_F, u_R, u, f_F, f_R) = SolveFEM(Kglob, u, f, restricted_DOF)
    if print_f_R
        println(f_R)
    end

    # Strains and stresses inside the element
    # print("Computing stress and strain..."); @time (strains_el, stresses_el) = StrainsStressesElem(v_nodes, v_mesh, E_mesh, pois_mesh, u)
    (strains_el, stresses_el) = StrainsStressesElem(v_nodes, v_mesh, E_mesh, pois_mesh, u)
    (voigtstrains_el, voigtstresses_el) = VoigtFormula(strains_el, stresses_el)
    # Strains and stresses calculated in each node
    (strains_nod, stresses_nod) = StrainsStressesNodes(nodal_connectivity, strains_el, stresses_el)

    return voigtstrains_el, voigtstresses_el, strains_nod, stresses_nod, strains_el, stresses_el

end
