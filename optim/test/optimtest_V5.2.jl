#= In this version, the optimization of the cube is performed in one single step,
without dividing the cube into different layers.
Version 5.2 add a limitation to the maximum stiffness value (Emax).
If the element's stiffness is > Emax, it is not updated. However, if E < Emax it
can be updated to a value > Emax, but in the next loop it will be removed from
the optimization.
=#

using Revise
using tfgfem
using optim
using Statistics
using Plots

cd("C:\\Users\\Hugo\\Desktop\\optim\\LOAD CASE $LC\\V5\\V5.2")

Emax = 10000.0

Eold2_mesh = copy(E_mesh)

compliance_collector = []
comp_count = 0

aux2collector = []
Hi2collector = []
Wi2collector = []
Ei2collector = []

elements_optim = []
for i = 1:elements
    push!(elements_optim, i)
end

loop2 = 0
conv_cont = 0
#while loop2 < 20
while conv_cont < 2
    global  Eold2_mesh, comp_count, loop2, conv_cont
    loop2 = loop2 + 1
    println("loop2 number ", loop2)

    num_elements_optim = length(elements_optim)

    (voigtstrains_el2, voigtstresses_el2, strains_nod2, stresses_nod2, strains_el2, stresses_el2) = FEAnalysis(v_nodes,
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
    for i = 1:num_elements_optim

        iel = elements_optim[i]

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
    #push!(workcollector2, work2)
    meanHi2 = mean(v_Hi2)
    varHi2 = var(v_Hi2)
    devHi2 = sqrt(varHi2)

    v_aux2 = []
    delete_vector = []
    for i = 1:num_elements_optim

        iel = elements_optim[i]

        aux2 = ((v_Hi2[i]-meanHi2))/devHi2 + 1
        if (aux2 < 0) || (Eold2_mesh[iel] >= Emax)
            println("negative aux or E>Emax at: ", i)
            push!(delete_vector, i)
        else
            push!(v_aux2, aux2)
            Enew2 = Eold2_mesh[iel]*aux2
            Eold2_mesh[iel] = Enew2
        end

    end
    push!(aux2collector, v_aux2)

    deleteat!(elements_optim, delete_vector)

    E2_nodes = NodalYoungModule(nodal_connectivity, Eold2_mesh)
    PrintParaview(v_nodes, v_mesh, strains_nod2, stresses_nod2, u, E2_nodes, "StiffEvolV52_$loop2.vtu")


    #loop2 = loop2 + 1
end #convergence loop2 (while)

(voigtstrains_el2, voigtstresses_el2, strains_nod2, stresses_nod2) =
FEAnalysis(v_nodes, v_mesh, Eold2_mesh, pois_mesh, u, f, restricted_DOF)

E2_nodes = NodalYoungModule(nodal_connectivity, Eold2_mesh)
PrintParaview(v_nodes, v_mesh, strains_nod2, stresses_nod2, u, E2_nodes, "OptimResultsV52($(loop2)).vtu")

plot(compliance_collector, title = "Compliance evolution. LC$LC (V5.2) (Emax=$Emax)", ylims = (0,210000))
#savefig("C:\\Users\\Hugo\\Desktop\\optim\\Compliance evolution. LC4 (V5.2) (Emax=$Emax).png")
