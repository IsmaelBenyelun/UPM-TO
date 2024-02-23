#=In this version, the optimization of the cube is performed in one single step,
without dividing the cube into different layers. The difference between V5.2 and
V5.1 is the following: in this version we strictly limit the maximum stiffness of
the elements, while in V5.2 elements where eliminated only if their stiffness was
over the maximum PRIOR to the update, not after.
=#

using Revise
using tfgfem
using optim
using Statistics
using Plots

cd("C:\\Users\\Hugo\\Desktop\\optim\\LOAD CASE $LC\\V5\\V5.3")

tick()
Emax = 20000.0

Eold_mesh = copy(E_mesh)

compliance_collector = []
comp_count = 0

aux_collector = []
Hi_collector = []
Wi_collector = []
Ei_collector = []

elements_optim = []
for i = 1:elements
    push!(elements_optim, i)
end

loop = 0
conv_cont = 0

while conv_cont < 2
    global  Eold_mesh, comp_count, loop, conv_cont
    loop = loop + 1
    println("loop number ", loop)

    num_elements_optim = length(elements_optim)

    (voigtstrains_el2, voigtstresses_el2, strains_nod2, stresses_nod2, strains_el2, stresses_el2) = FEAnalysis(v_nodes,
    v_mesh, Eold_mesh, pois_mesh, u, f, restricted_DOF)
    (K, B) = Kglobal(dof_node, v_nodes, v_mesh, Eold_mesh, pois_mesh)

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
    if conv_cont == 2
        loop = loop - 1
        break
    end

    v_Hi = []
    v_Wi = []
    v_Ei = []
    for i = 1:num_elements_optim

        iel = elements_optim[i]

        push!(v_Ei, Eold_mesh[iel])

        Dr = DrMatrix(pois)
        strains = voigtstrains_el2[iel]
        elem_coord = CoordsElement(v_mesh, v_nodes, iel)
        vol_iel = VolumeTetra(elem_coord)

        Hi = ConstantH(Dr, strains, vol_iel)
        if Hi < 0
            println("negative H in element: ", iel)
        end
        push!(v_Hi, Hi)

        Wi = Eold_mesh[iel] * Hi
        push!(v_Wi, Wi)

    end
    push!(Hi_collector, v_Hi)
    push!(Wi_collector, v_Wi)
    push!(Ei_collector, v_Ei)

    work = sum(v_Wi)
    #push!(workcollector2, work)
    meanHi = mean(v_Hi)
    varHi = var(v_Hi)
    devHi = sqrt(varHi)

    v_aux = []
    delete_vector = []
    for i = 1:num_elements_optim

        iel = elements_optim[i]

        aux = ((v_Hi[i]-meanHi))/devHi + 1
        if (aux < 0) || (Eold_mesh[iel] >= Emax)
            println("negative aux or E>Emax at: ", i)
            push!(delete_vector, i)
        else
            push!(v_aux, aux)
            Enew = Eold_mesh[iel]*aux
            if Enew >= Emax
                Eold_mesh[iel] = Emax
            else
                Eold_mesh[iel] = Enew
            end

        end

    end
    push!(aux_collector, v_aux)

    deleteat!(elements_optim, delete_vector)

    E_nodes = NodalYoungModule(nodal_connectivity, Eold_mesh)
    PrintParaview(v_nodes, v_mesh, strains_nod2, stresses_nod2, u, E_nodes, "StiffEvolV53_$loop.vtu")

    #loop = loop + 1
end #convergence loop (while)

(voigtstrains_el2, voigtstresses_el2, strains_nod2, stresses_nod2) =
FEAnalysis(v_nodes, v_mesh, Eold_mesh, pois_mesh, u, f, restricted_DOF)

E_nodes = NodalYoungModule(nodal_connectivity, Eold_mesh)
PrintParaview(v_nodes, v_mesh, strains_nod2, stresses_nod2, u, E_nodes, "OptimResultsV53($(loop)).vtu")

tock()

plot(compliance_collector, title = "", label = "",
ylabel = "Compliance", xlabel = "Optimization Steps", ylims = (0,3000), xlims = (0,20))
#savefig("C:\\Users\\Hugo\\Desktop\\optim\\LOAD CASE $LC\\optim process\\Compliance evolution. LC$LC (Emax=$Emax).png")
