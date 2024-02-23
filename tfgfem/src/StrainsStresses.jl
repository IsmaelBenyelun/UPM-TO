"""
    (strains_el, stresses_el) = StrainsStressesElem(v_nodes, v_mesh, E_mesh, pois_mesh, u)

Compute the element's strains and stresses
u: displacement calculated in previous steps

"""

function StrainsStressesElem(v_nodes, v_mesh, E_mesh, pois_mesh, u)

    elements = length(v_mesh)

    strains_el=[]
    stresses_el=[]

    strains_el_eig=[]
    stresses_el_eig=[]
    displacements_el=[]
    coord_nodes_el=[]

    for i in range(1,stop=elements)

        el_loc=i
        E_el=E_mesh[i]
        pois_el = pois_mesh[i]

        # --- Mechanical model
        aux_vec_coord=[]
        disp_loc=[]
        #For each element i, j2 counter chooses one of its nodes
        for j2 in range(1,stop=4)
            naux=v_mesh[el_loc][j2]
            cont = 0
            # k2 counter selects each of the DOF from node j2
            for k2 = 2:-1:0
                cont = cont + 1
                push!(aux_vec_coord, v_nodes[naux][cont])
                push!(disp_loc, u[naux*3-k2])
            end
        end

        push!(displacements_el, disp_loc)
        push!(coord_nodes_el, aux_vec_coord)

        # VOLVEMOS A LLAMAR A LA FUNCION AUXILIAR DE LA RIGIDEZ LOCAL
        (Kloc, B, voltet, D, straM)=KlocTetra(aux_vec_coord, E_el, pois_el)

        stra_mat=B
        stress_mat=straM

        local_disp=disp_loc

        matloc=stra_mat
        vec_strain_loc=matloc*local_disp

        matloc=stress_mat
        vec_stress_loc=matloc*local_disp  #formulation: (εx, εy, εz, γxy, γyz, γxz)

        stress=[vec_stress_loc[1] vec_stress_loc[4] vec_stress_loc[6];
        vec_stress_loc[4] vec_stress_loc[2] vec_stress_loc[5];
        vec_stress_loc[6] vec_stress_loc[5] vec_stress_loc[3]]

#=
        stress=[vec_stress_loc[1] vec_stress_loc[4] vec_stress_loc[5];
        vec_stress_loc[4] vec_stress_loc[2] vec_stress_loc[6];
        vec_stress_loc[5] vec_stress_loc[6] vec_stress_loc[3]]
=#
        eigvals_stress=eigvals(stress)
        eig1 = eigvals_stress[1]
        eig2 = eigvals_stress[2]
        eig3 = eigvals_stress[3]
        push!(stresses_el_eig, [eig1,eig2,eig3])
        push!(stresses_el, stress)
        #stresses_el.append([vec_stress_loc[0],vec_stress_loc[1],vec_stress_loc[2]])

        strain =[vec_strain_loc[1] vec_strain_loc[4] vec_strain_loc[6];
        vec_strain_loc[4] vec_strain_loc[2] vec_strain_loc[5];
        vec_strain_loc[6] vec_strain_loc[5] vec_strain_loc[3]]

#=
        strain =[vec_strain_loc[1] vec_strain_loc[4] vec_strain_loc[5];
        vec_strain_loc[4] vec_strain_loc[2] vec_strain_loc[6];
        vec_strain_loc[5] vec_strain_loc[6] vec_strain_loc[3]]
=#

        eigvals_strain=eigvals(strain)
        eig1 = eigvals_strain[1]
        eig2 = eigvals_strain[2]
        eig3 = eigvals_strain[3]
        push!(strains_el_eig, [eig1,eig2,eig3])
        push!(strains_el, strain)
    end

    return (strains_el, stresses_el, strains_el_eig, stresses_el_eig) #Autovalores

end


"""
    (strains_nod, stresses_nod) = StrainsStressesNodes(nodal_connectivity, strains_el, stresses_el)

Compute the equivalent strains and stresses in the nodes, using the connectivity.
This must be done for further printing in programs such as ParaView

"""

function StrainsStressesNodes(nodal_connectivity, strains_el, stresses_el)

    nodes = length(nodal_connectivity)
    stresses_nod = []

    # recorremos cada nodo
    for i in range(1, stop=nodes)
        sigmax = 0
        sigmay = 0
        sigmaz = 0
        sigmaxy = 0
        sigmayz = 0
        sigmaxz = 0
        #recorremos cada elemento asociado a ese nodo
        for j in range(1, stop=length(nodal_connectivity[i]))
            sigmax += (stresses_el[nodal_connectivity[i][j]][1,1]) / (length(nodal_connectivity[i]))
            sigmay += (stresses_el[nodal_connectivity[i][j]][2,2]) / (length(nodal_connectivity[i]))
            sigmaz += (stresses_el[nodal_connectivity[i][j]][3,3]) / (length(nodal_connectivity[i]))
            sigmaxy += (stresses_el[nodal_connectivity[i][j]][1,2]) / (length(nodal_connectivity[i]))
            sigmayz += (stresses_el[nodal_connectivity[i][j]][2,3]) / (length(nodal_connectivity[i]))
            sigmaxz += (stresses_el[nodal_connectivity[i][j]][1,3]) / (length(nodal_connectivity[i]))
        end
        push!(stresses_nod, [sigmax, sigmay, sigmaz, sigmaxy, sigmayz, sigmaxz])
    end

    strains_nod = []
    # recorremos cada nodo
    for i in range(1, stop=nodes)
        epsx = 0
        epsy = 0
        epsz = 0
        epsxy = 0
        epsyz = 0
        epsxz = 0
        #recorremos cada elemento asociado a ese nodo
        for j in range(1, stop=length(nodal_connectivity[i]))
            epsx += (strains_el[nodal_connectivity[i][j]][1,1]) / (length(nodal_connectivity[i]))
            epsy += (strains_el[nodal_connectivity[i][j]][2,2]) / (length(nodal_connectivity[i]))
            epsz += (strains_el[nodal_connectivity[i][j]][3,3]) / (length(nodal_connectivity[i]))
            epsxy += (strains_el[nodal_connectivity[i][j]][1,2]) / (length(nodal_connectivity[i]))
            epsyz += (strains_el[nodal_connectivity[i][j]][2,3]) / (length(nodal_connectivity[i]))
            epsxz += (strains_el[nodal_connectivity[i][j]][1,3]) / (length(nodal_connectivity[i]))


        end
        push!(strains_nod, [epsx, epsy, epsz, epsxy, epsyz, epsxz])
    end

    return (strains_nod, stresses_nod)

end


"""
    (voigt_strains, voigt_stresses) = VoigtFormula(strains_el, stresses_el)

Rewrites the element's strains and stresses in the Voigt formulation.
The order of the strains terms is the following: (εx, εy, εz, γxy, γyz, γxz)
And the same for the stresses.
"""
function VoigtFormula(strains_el, stresses_el)

    elements = length(strains_el)

    voigt_strains = []
    voigt_stresses = []

    for i = 1:elements

        # (εx, εy, εz, γxy, γyz, γxz)
        push!(voigt_strains, [strains_el[i][1,1], strains_el[i][2,2], strains_el[i][3,3],
        strains_el[i][1,2], strains_el[i][2,3], strains_el[i][1,3]])

        push!(voigt_stresses, [stresses_el[i][1,1], stresses_el[i][2,2], stresses_el[i][3,3],
        stresses_el[i][1,2], stresses_el[i][2,3], stresses_el[i][1,3]])

    end

    return voigt_strains, voigt_stresses

end
