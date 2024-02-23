#In this version, V2, I corrected the issue with the absolute value but I began to
#to obtain negative values for the stiffness.
using Revise
using tfgfem
using optim
using Statistics
using Plots

Eold2_mesh = copy(E_mesh)
Enew2_mesh = copy(E_mesh)
# Layer we are analyzing
#for lay = 10:-1:9
aux2collector = []
Hi2collector = []
Wi2collector = []
Ei2collector = []
lay = 10
    println("analyzing layer ", lay)
    num_elements_layer = length(layers[lay])  #number of elements in the layer
    elements_layer = layers[lay]   #collection of elements in the layer

    workcollector2 = []
    loop2 = 1
    while loop2 <= 20
        global  Eold2_mesh, strains_nod2, stresses_nod2, loop2
        println("loop2 number ", loop2)

        (voigtstrains_el2, voigtstresses_el2, strains_nod2, stresses_nod2) = FEAnalysis(v_nodes,
        v_mesh, Eold2_mesh, pois_mesh, u, f, restricted_DOF)

        v_Hi2 = []
        v_Wi2 = []
        v_Ei2 = []
        for i = 1:num_elements_layer

            iel = layers[lay][i]

            push!(v_Ei2, Eold2_mesh[iel])

            Dr = DrMatrix(pois)
            strains = voigtstrains_el2[iel]
            elem_coord = CoordsElement(v_mesh, v_nodes, iel)
            vol_iel = VolumeTetra(elem_coord)

            Hi2 = ConstantH(Dr, strains, vol_iel)
            if Hi2 < 0
                println("negative H in element: ", iel)
            end
            push!(v_Hi2, Hi2)

            Wi2 = Eold2_mesh[iel] * Hi2
            push!(v_Wi2, Wi2)

        end
        push!(Hi2collector, v_Hi2)
        push!(Wi2collector, v_Wi2)
        push!(Ei2collector, v_Ei2)

        work2 = sum(v_Wi2)
        push!(workcollector2, work2)
        meanHi2 = mean(v_Hi2)
        varHi2 = var(v_Hi2)
        devHi2 = sqrt(varHi2)

        v_aux2 = []
        for i = 1:num_elements_layer

            iel = layers[lay][i]
            aux2 = ((v_Hi2[i]-meanHi2))/devHi2 + 1
            push!(v_aux2, aux2)
            Enew2 = Eold2_mesh[iel]*aux2
            Enew2_mesh[iel] = Enew2
        end
        push!(aux2collector, v_aux2)

        #E2_nodes = NodalYoungModule(nodal_connectivity, Eold2_mesh)
        #PrintParaview(v_nodes, v_mesh, strains_nod2, stresses_nod2, u, E2_nodes, "StiffnessEvolution_layer$(lay)_$loop2.vtu")

        Eold2_mesh = copy(Enew2_mesh)

        loop2 = loop2 + 1
    end #convergence loop2 (while)

    (voigtstrains_el2, voigtstresses_el2, strains_nod2, stresses_nod2) =
    FEAnalysis(v_nodes, v_mesh, Enew2_mesh, pois_mesh, u, f, restricted_DOF)
    E2_nodes = NodalYoungModule(nodal_connectivity, Eold2_mesh)
    PrintParaview(v_nodes, v_mesh, strains_nod2, stresses_nod2, u, E2_nodes, "OptimResults2_layer$(lay)($(loop2-1)).vtu")

#end #layer loop

E2_mesh = copy(Enew2_mesh)
(voigtstrains_el2, voigtstresses_el2, strains_nod2, stresses_nod2) =
FEAnalysis(v_nodes, v_mesh, E2_mesh, pois_mesh, u, f, restricted_DOF)

nodal_connectivity = NodalConnectivity(v_nodes, v_mesh)
E_nodes = NodalYoungModule(nodal_connectivity, E2_mesh)
PrintParaview(v_nodes, v_mesh, strains_nod2, stresses_nod2, u, E_nodes, "mesh_results2_$(lay).vtu")
