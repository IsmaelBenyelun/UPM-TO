#In this version 4, I am going to try to fix the optimization for the
#next layers. I am going to implent code that so that the following layer
#starts with stiffness values from the previous layer.

using Revise
using tfgfem
using optim
using Statistics
using Plots

cd("C:\\Users\\Hugo\\Desktop\\optim\\V4")

Eold2_mesh = copy(E_mesh)

# Layer we are analyzing
compliance_collector = []
comp_count = 0
for lay = 10:-1:1
#lay = 9
    aux2collector = []
    Hi2collector = []
    Wi2collector = []
    Ei2collector = []
    println("analyzing layer ", lay)

    elements_layer = copy(layers[lay])   #elements that will be optimize
    if lay != 10
        for i = 1:length(layers[lay])
            iel = layers[lay][i]
            Eold2_mesh[iel] = copy(Emeshlayer_collector[lay+1][i])
        end
    end

    workcollector2 = []
    loop2 = 1
    conv_cont = 0
    while conv_cont < 2
        global  Eold2_mesh, comp_count#, loop2, conv_cont
        println("loop2 number ", loop2)

        num_elements_layer = length(elements_layer)

        (voigtstrains_el2, voigtstresses_el2, strains_nod2, stresses_nod2) = FEAnalysis(v_nodes,
        v_mesh, Eold2_mesh, pois_mesh, u, f, restricted_DOF)
        (K, B) = Kglobal(dof_node, v_nodes, v_mesh, Eold2_mesh, pois_mesh)

        uT = transpose(u)
        compliance = 0.5 * uT * K * u
        push!(compliance_collector, compliance)
        comp_count = comp_count + 1
        #Convergence criteria
        if comp_count > 1
            diff = abs((compliance_collector[comp_count]-compliance_collector[comp_count-1]))/compliance_collector[comp_count-1]
            if diff < 0.005
                conv_cont = conv_cont + 1
            else
                conv_cont = 0
            end
        end

        v_Hi2 = []
        v_Wi2 = []
        v_Ei2 = []
        for i = 1:num_elements_layer

            iel = elements_layer[i]

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
        delete_vector = []
        for i = 1:num_elements_layer

            iel = elements_layer[i]

            aux2 = ((v_Hi2[i]-meanHi2))/devHi2 + 1
            if aux2 < 0
                println("negative aux at: ", i)
                push!(delete_vector, i)
            else
                push!(v_aux2, aux2)
                Eold2_mesh[iel] = Eold2_mesh[iel] * aux2
            end

        end
        push!(aux2collector, v_aux2)

        deleteat!(elements_layer, delete_vector)

        E2_nodes = NodalYoungModule(nodal_connectivity, Eold2_mesh)
        PrintParaview(v_nodes, v_mesh, strains_nod2, stresses_nod2, u, E2_nodes, "StiffEvolV4_layer$(lay)_$loop2.vtu")

        loop2 = loop2 + 1
    end #convergence loop2 (while)

    for i = 1:length(Emeshlayer_collector[lay])
        iel = layers[lay][i]
        Emeshlayer_collector[lay][i] = copy(Eold2_mesh[iel])
    end

    (voigtstrains_el2, voigtstresses_el2, strains_nod2, stresses_nod2) =
    FEAnalysis(v_nodes, v_mesh, Eold2_mesh, pois_mesh, u, f, restricted_DOF)

    E2_nodes = NodalYoungModule(nodal_connectivity, Eold2_mesh)
    PrintParaview(v_nodes, v_mesh, strains_nod2, stresses_nod2, u, E2_nodes, "OptimResultsV4_layer$(lay)($(loop2-1)).vtu")

end #layers loop
