module Optimfunction

using optim
using tfgfem

using Statistics
using LinearAlgebra
using DelimitedFiles

export optimfunction

function optimfunction(;E_max, k, tol, E_mesh, elements, v_nodes, v_mesh, ν_mesh, u, f, restricted_DOF,
dof_node, nodal_connectivity, ν, vol, save_csv=false)

    # Storage
    compliance_collector = []
    aux_collector = []
    Hi_collector = []
    Wi_collector = []
    Ei_collector = []
    elements_optim = [i for i in 1:elements]

    # Young modulus and energy mesh
    Eold_mesh = copy(E_mesh)
    α_mesh = zeros(elements)
    W_mesh = zeros(elements)

    # Convergence counter
    loop = 0
    conv_count = 0
    comp_count = 0

    while conv_count < 2
        """INITIAL FINITE ELEMENT ANALYSIS"""
        global voigtstrains_el2, v_Wi
        if loop == 0
            voigtstrains_el2, voigtstresses_el2, strains_nod2, stresses_nod2, strains_el2, stresses_el2 = FEAnalysis(v_nodes,
            v_mesh, Eold_mesh, ν_mesh, u, f, restricted_DOF)

            # K, B = Kglobal(dof_node, v_nodes, v_mesh, Eold_mesh, ν_mesh)
            # uT = transpose(u)
            # compliance = 0.5 * uT * K * u
            compliance = 0.5 * f ⋅ u
            push!(compliance_collector, compliance)
            comp_count = comp_count + 1

            E_nodes = NodalYoungModule(nodal_connectivity, Eold_mesh)
            α_nodes = NodalYoungModule(nodal_connectivity, α_mesh)
            W_nodes = NodalYoungModule(nodal_connectivity, W_mesh)

            PrintParaview(v_nodes, v_mesh, strains_nod2, stresses_nod2, u, E_nodes, α_nodes, α_mesh, W_nodes, "./paraview-output/StiffEvolV53_$loop.vtu", verbose=false)
            println("Initial point. Compliance: $compliance")
        end

        loop += 1
        #println("loop number ", loop)
        # num_elements_optim = length(elements_optim) # pythonic update

        """STRAIN PARAMETERS Hi COMPUTATION"""
        v_Hi = []
        v_Wi = []
        v_Ei = []

        for iel in elements_optim # pythonic update
        # for i = 1:num_elements_optim
            # iel = elements_optim[i]

            push!(v_Ei, Eold_mesh[iel])

            Dr = DrMatrix(ν)
            strains = voigtstrains_el2[iel]
            vol_iel = vol[iel]

            Hi = ConstantH(Dr, strains, vol_iel)
            push!(v_Hi, Hi)

            Wi = Eold_mesh[iel] * Hi
            W_mesh[iel] = Wi
            push!(v_Wi, Wi)

        end

        push!(Hi_collector, v_Hi)
        push!(Wi_collector, v_Wi)
        push!(Ei_collector, v_Ei)

        """YOUNG MODULUS Ei UPDATE"""
        meanHi = mean(v_Hi)
        varHi = var(v_Hi)
        devHi = sqrt(varHi)

        v_aux = []
        delete_vector = []

        for (i, iel) in enumerate(elements_optim)

            α = ((v_Hi[i]-meanHi))/(k * devHi)

            aux = 1 + α
            α_mesh[iel] = α
            if (aux < 0) || (Eold_mesh[iel] >= E_max) # || is 'or'
                push!(delete_vector, i) # Apply deleteat! elements optim
            else
                push!(v_aux, aux)
                Enew = Eold_mesh[iel] * aux

                if Enew >= E_max
                    Eold_mesh[iel] = E_max

                else
                    Eold_mesh[iel] = Enew
                end
            end
        end

        push!(aux_collector, v_aux)
        deleteat!(elements_optim, delete_vector) # Important!!! Does not delete the element but the POSITION
        # filter!(e -> e ∉ delete_vector, elements_optim) # Deletes the element!
        # append!(delete_vector_collector, delete_vector)

        """CONVERGENCE ANALYSIS"""
        voigtstrains_el2, voigtstresses_el2, strains_nod2, stresses_nod2, strains_el2, stresses_el2 = FEAnalysis(v_nodes,
        v_mesh, Eold_mesh, ν_mesh, u, f, restricted_DOF)
        E_nodes = NodalYoungModule(nodal_connectivity, Eold_mesh)
        α_nodes = NodalYoungModule(nodal_connectivity, α_mesh)
        W_nodes = NodalYoungModule(nodal_connectivity, W_mesh)

        PrintParaview(v_nodes, v_mesh, strains_nod2, stresses_nod2, u, E_nodes, α_nodes, α_mesh, W_nodes, "./paraview-output/StiffEvolV53_$loop.vtu", verbose=false)

        compliance = 0.5 * f ⋅ u

        push!(compliance_collector, compliance)
        comp_count = comp_count + 1
        # Verbosity
        println("Iteration $loop. Compliance: $compliance")
        if comp_count > 1
            diff = abs((compliance_collector[comp_count]-compliance_collector[comp_count-1]))/compliance_collector[comp_count-1]
            if diff < tol
                conv_count = conv_count + 1
            else
                conv_count = 0
            end
        end
    end # convergence loop (while)

    #println("loop number ", loop)

    # Printing (commented in order to optimize code runtime)
    (voigtstrains_el2, voigtstresses_el2, strains_nod2, stresses_nod2) =
    FEAnalysis(v_nodes, v_mesh, Eold_mesh, ν_mesh, u, f, restricted_DOF)#, print_f_R=true)

    # Final volume
    ρ = Eold_mesh ./ E_max
    volfrac = transpose(ρ) * vol / sum(vol)

    filter = true ###

    # println(elements_optim)
    if filter
        all_nodes_optim = reduce(vcat, v_mesh[elements_optim])
        nodes_optim = unique(sort(all_nodes_optim))
        # nodal_connectivity = MeshFilter(nodal_connectivity, elements_optim)
    end
    E_nodes = NodalYoungModule(nodal_connectivity, Eold_mesh)
    α_nodes = NodalYoungModule(nodal_connectivity, α_mesh)
    W_nodes = NodalYoungModule(nodal_connectivity, W_mesh)

    PrintParaview(v_nodes, v_mesh, strains_nod2, stresses_nod2, u, E_nodes, α_nodes, α_mesh, W_nodes, "./paraview-output/OptimResultsV53_last-v=$volfrac.vtu", verbose=true)

    if save_csv
        writedlm("./paraview-output/complianceV53.csv", compliance_collector, ",")
        writedlm("./paraview-output/E_mesh.csv", Eold_mesh, ",")
        writedlm("./paraview-output/W_mesh.csv", W_mesh, ",")
        writedlm("./paraview-output/W_nodes.csv", W_nodes, ",")
        writedlm("./paraview-output/elements.csv", elements_optim, ",")
        writedlm("./paraview-output/nodes.csv", nodes_optim, ",")
    end

    return Eold_mesh, E_nodes, volfrac
end

end # MODULE
