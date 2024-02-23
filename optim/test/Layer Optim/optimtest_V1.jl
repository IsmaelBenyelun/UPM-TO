#This is the first version of the optimizing function. I made the mistake of using
#an absolute value where it didn't correspond
using Revise
using tfgfem
using optim
using Statistics
using Plots

Eold_mesh = copy(E_mesh)
Enew_mesh = copy(E_mesh)
# Layer we are analyzing
#for lay = 10:-1:1
auxcollector = []
Hicollector = []
Wicollector = []
Eicollector = []
lay = 10
    println("analyzing layer ", lay)
    elements_layer = length(layers[lay])

    workcollector = []
    loop = 1

    while loop <= 30
    #while abs(mean_constr_min) >= 0.001
        global  Eold_mesh, loop
        println("loop number ", loop)

        (voigtstrains_el, voigtstresses_el, strains_nod, stresses_nod) = FEAnalysis(v_nodes,
        v_mesh, Eold_mesh, pois_mesh, u, f, restricted_DOF)

        v_Hi = []
        v_Wi = []
        v_Ei = []
        for i = 1:elements_layer

            iel = layers[lay][i]

            push!(v_Ei, Eold_mesh[iel])

            Dr = DrMatrix(pois)
            strains = voigtstrains_el[iel]
            elem_coord = CoordsElement(v_mesh, v_nodes, iel)
            vol_iel = VolumeTetra(elem_coord)

            H = ConstantH(Dr, strains, vol_iel)
            if H < 0
                println("negative H in element: ", iel)
            end
            push!(v_Hi, H)

            Wi = Eold_mesh[iel] * H
            push!(v_Wi, Wi)


        end
        push!(Hicollector, v_Hi)
        push!(Wicollector, v_Wi)
        push!(Eicollector, v_Ei)

        work = sum(v_Wi)
        push!(workcollector, work)

        meanH = mean(v_Hi)
        varH = var(v_Hi)
        devH = sqrt(varH)

        v_aux = []
        for i = 1:elements_layer

            iel = layers[lay][i]
            aux = (abs(v_Hi[i]-meanH))/devH + 1
            push!(v_aux, aux)

            Enew = Eold_mesh[iel]*aux
            Enew_mesh[iel] = Enew
        end
        push!(auxcollector, v_aux)

        #E_nodes = NodalYoungModule(nodal_connectivity, Eold_mesh)
        #PrintParaview(v_nodes, v_mesh, strains_nod, stresses_nod, u, E_nodes, "StiffnessEvolution$loop.vtu")

        Eold_mesh = copy(Enew_mesh)

        loop = loop + 1
    end #convergence loop (while)

#end #layer loop

E_mesh = copy(Enew_mesh)
(voigtstrains_el, voigtstresses_el, strains_nod, stresses_nod) = FEAnalysis(v_nodes, v_mesh, E_mesh, pois_mesh, u, f, restricted_DOF)

nodal_connectivity = NodalConnectivity(v_nodes, v_mesh)
E_nodes = NodalYoungModule(nodal_connectivity, E_mesh)
PrintParaview(v_nodes, v_mesh, strains_nod, stresses_nod, u, E_nodes, "mesh_results_$lay.vtu")
