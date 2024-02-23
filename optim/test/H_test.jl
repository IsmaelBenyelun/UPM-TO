using Revise
using tfgfem
using optim
using Statistics


oldE_mesh = E_mesh
v_Hplot = zeros(length(E_mesh))
# Layer we are analyzing
#for lay = 10:-1:1
lay = 1
    println("analyzing layer ", lay)
    elements_layer = length(layers[lay])

    i = 1

    #while i <= 6
    #while abs(mean_constr_min) >= 0.001
        #global  oldE_mesh, mean_constr_min, v_H, i
        println("loop number ", i)

        (voigtstrains_el, voigtstresses_el, strains_nod, stresses_nod) = FEAnalysis(v_nodes,
        v_mesh, oldE_mesh, pois_mesh, u, f, restricted_DOF)

        v_H = []
        for i = 1:elements_layer

            iel = layers[lay][i]
            Dr = DrMatrix(pois)
            strains = voigtstrains_el[iel]
            elem_coord = CoordsElement(v_mesh, v_nodes, iel)
            vol_iel = VolumeTetra(elem_coord)

            H = ConstantH(Dr, strains, vol_iel)
            if H < 0
                println("negative H in element: ", iel)
            end
            push!(v_H, H)
            v_Hplot[iel] = H


        end

        meanH = mean(v_H)
        varH = var(v_H)
        devH = sqrt(varH)

        i = i + 1
    #end #convergence loop (while)

#end #layer loop


nodal_connectivity = NodalConnectivity(v_nodes, v_mesh)
H_nodes = NodalYoungModule(nodal_connectivity, v_Hplot)
PrintParaview(v_nodes, strains_nod, stresses_nod, u, H_nodes, v_mesh)

#v_Hplotnodes = NodalYoungModule(nodal_connectivity, v_Hplot)
#PrintParaview(v_nodes, strains_nod, stresses_nod, u, v_Hplotnodes,v_mesh)
