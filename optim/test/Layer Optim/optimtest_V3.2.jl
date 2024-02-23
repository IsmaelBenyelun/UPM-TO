#=In V3 the negative values for stiffess have been corrected by deleting
from the optimization those elements that reach negative values.
However, there are still problems for the layers following no.10
To try to solve this problem, a maximum stiffness is going to be imposed.
As a result, elements with E>Emax won't be updated.
=#
using Revise
using tfgfem
using optim
using Statistics
using Plots

cd("C:\\Users\\Hugo\\Desktop\\optim\\V3\\V3.2")

Emax = 6000.0

Eold2_mesh = copy(E_mesh)
# Layer we are analyzing
compliance_collector = []
comp_count = 0
for lay = 10:-1:1
aux2collector = []
Hi2collector = []
Wi2collector = []
Ei2collector = []
#lay = 7
    println("analyzing layer ", lay)
    #num_elements_layer = length(layers[lay])  #number of elements in the layer
    elements_layer = copy(layers[lay])   #collection of elements remaining to optimize in layer

    workcollector2 = []
    #compliance_collector = []
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
            if (aux2 < 0) || (Eold2_mesh[iel] >= Emax)
                println("negative aux or E>Emax at: ", i)
                push!(delete_vector, i)
            else
                push!(v_aux2, aux2)
                Eold2_mesh[iel]  = Eold2_mesh[iel]*aux2
            end

        end
        push!(aux2collector, v_aux2)

        deleteat!(elements_layer, delete_vector)

        E2_nodes = NodalYoungModule(nodal_connectivity, Eold2_mesh)
        PrintParaview(v_nodes, v_mesh, strains_nod2, stresses_nod2, u, E2_nodes, "StiffEvolV32_layer$(lay)_$loop2.vtu")

        Eold2_mesh = copy(Eold2_mesh)

        loop2 = loop2 + 1
    end #convergence loop2 (while)

    (voigtstrains_el2, voigtstresses_el2, strains_nod2, stresses_nod2) =
    FEAnalysis(v_nodes, v_mesh, Eold2_mesh, pois_mesh, u, f, restricted_DOF)

    E2_nodes = NodalYoungModule(nodal_connectivity, Eold2_mesh)
    PrintParaview(v_nodes, v_mesh, strains_nod2, stresses_nod2, u, E2_nodes, "OptimResultsV32_layer$(lay)($(loop2-1)).vtu")

end #layers loop
