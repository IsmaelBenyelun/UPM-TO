"""Script to run and understand Generative Design scripts developed by
Saucedo et al. in "Generative design methodology using an alternative
topology optimization based on stiffness modification"
"""

push!(LOAD_PATH, pwd())
using tfgfem
using Optimfunction
using DelimitedFiles

# GLOBAL VARIABLES - Properties
E0 = 200.0e3 # MPa
ν = 0.33
E_MAX = 500.0e3 # MPa

# GLOBAL VARIABLES
DOF_NODE = 3

SIZE_X = 100. # mm
SIZE_Y = 100. # mm
SIZE_Z = 100. # mm

GAP_X =  12.5 # 5. mm
GAP_Y = 12.5 # 5. mm
GAP_Z = 12.5 # 5. mm

# GLOBAL VARIABLES - DIRICHLET BC
OUTER_RADIUS_BC = 35.0 # mm
INNER_RADIUS_BC = 25.0 # mm

RADIUS_BC = 20.0 # mm
FORCE_MAGNITUDE = 1.0e4 # N

# UPDATE PARAMETER
K = 1

# TOLERANCE
TOL = 1e-2

# LOAD CASE (Un)comment the one to run
LOAD_CASE = "TUBE"
# LOAD_CASE = "SHEAR"
# LOAD_CASE = "CHAIR"
# LOAD_CASE = "TORSION"

function build_mesh(size_x=SIZE_X, size_y=SIZE_Y, size_z=SIZE_Z, gap_x=GAP_X, gap_y=GAP_Y, gap_z=GAP_Z, jump_x=0., dof_node=DOF_NODE; print_paraview=false)
    # Create array of mesh size and discretization
    v_size = [size_x, size_y, size_z]
    v_gap = [gap_x, gap_y, gap_z]
    # Obtaining vector por nodes, mesh and mechanical properties
    (nodes, elements, v_nodes, v_mesh, v_mechanical) = MeshBuild(v_size, v_gap, E0, ν)
    (E_mesh, ν_mesh) = MechanicalMesh(v_mechanical)
    mesh_distribution = MeshDistribution(v_size, v_gap, jump_x)
    # Number of degrees of freedom and nodal connectivity
    num_DOF = nodes*dof_node
    nodal_connectivity = NodalConnectivity(v_nodes, v_mesh)
    # Create vtu file with the mesh
    if print_paraview
        if !isdir("paraview-output")
            mkdir("paraview-output")
        end
        PrintParaview(v_nodes, v_mesh, "./paraview-output/mesh.vtu")
    end
    vol = [VolumeTetra(CoordsElement(v_mesh, v_nodes, i)) for i in 1:elements]
    return nodes, elements, v_nodes, v_mesh, E_mesh, ν_mesh, mesh_distribution, num_DOF, nodal_connectivity, vol
end

function set_boundary_conditions(v_nodes, mesh_distribution, num_DOF)
    # Set Dirichlet Boundary Conditions (i.e. displacements)
    restricted_DOF, u = _set_dirichlet_boundary_conditions(v_nodes, mesh_distribution, num_DOF, load_case=LOAD_CASE)
    # Set Neumann Boundary Conditions (i.e. *nodal* forces)
    force_DOF, f = _set_neumann_boundary_conditions(v_nodes, mesh_distribution, num_DOF, load_case=LOAD_CASE)
    return restricted_DOF, u, force_DOF, f
end

function run_optimization(elements, v_nodes, v_mesh, E_mesh, ν_mesh, u, f, restricted_DOF, nodal_connectivity, vol, ν=ν, dof_node=DOF_NODE, E_max=E_MAX, tol=TOL)
    return optimfunction(
        E_max=E_max, k=K, tol=tol, dof_node=dof_node,
        elements=elements, v_nodes=v_nodes, v_mesh=v_mesh, nodal_connectivity=nodal_connectivity,
        E_mesh=E_mesh, ν_mesh=ν_mesh, ν=ν,
        u=u, f=f, restricted_DOF=restricted_DOF, vol=vol,
        save_csv=true
    )
end

function _set_dirichlet_boundary_conditions(v_nodes, mesh_distribution, num_DOF, outer_radius_BC=OUTER_RADIUS_BC, inner_radius_BC=INNER_RADIUS_BC, radius_BC=RADIUS_BC; load_case=LOAD_CASE)
    """"SETTING BOUNDARY CONDITIONS: to see what the optimization is going to be"""
    conditionBC = [1, 1, 1] # All the displacements are restricted (applies on the 3 directions)
    BC = [0.0, 0.0, 0.0] # Value of the restrictions
    BCMat = [] # Initialize BCMat
    if issubset("TUBE", load_case) || load_case == "TORSION"
        # BOUNDARY CONDITIONS. EMBEDDED CIRCLE AT THE BOTTOM OF THE CUBE

        # Selecting affected nodes (at the center)
        x_bc = mesh_distribution[1][Int((length(mesh_distribution[1]) + 1)/2)]
        y_bc = mesh_distribution[2][Int((length(mesh_distribution[1]) + 1)/2)]
        z_bc = mesh_distribution[3][1]

        bc_coords = [x_bc, y_bc, z_bc]

        # CircumferenceZ will be the vector containing the nodes affected by the future BC.
        nodes_BC = CircunferenceZ(v_nodes, bc_coords, inner_radius_BC, outer_radius_BC)

    elseif load_case == "SHEAR"
        # Selecting affected nodes
        x_bc = mesh_distribution[1][1]
        y_bc = mesh_distribution[2][Int((length(mesh_distribution[1]) + 1)/2)]
        z_bc = mesh_distribution[3][Int((length(mesh_distribution[1]) + 1)/2)]
        bc_coords = [x_bc, y_bc, z_bc]

        radiusBC = 20.0
        # CircleZ will be the vector containing the nodes affected by the future BC.
        indicator = [1, 0, 0]
        nodes_BC = Circle(v_nodes, indicator, bc_coords, radiusBC)

    elseif load_case == "CHAIR" ## TODO: Check the center of the nodes!!!!!
        radius_disk = 5.
        nodes_BC = []
        #Node 1
        x_n1 = mesh_distribution[1][5]
        y_n1 = mesh_distribution[2][5]
        z_n1 = mesh_distribution[3][1]
        n1_coords = [x_n1, y_n1, z_n1]
        # n1 = NodeSelector(v_nodes, n1_coords)
        n1 = CircleZ(v_nodes, n1_coords, radius_disk)
        append!(nodes_BC, n1)

        #Node 2
        x_n2 = mesh_distribution[1][end-4]
        y_n2 = mesh_distribution[2][5]
        z_n2 = mesh_distribution[3][1]
        n2_coords = [x_n2, y_n2, z_n2]
        # n2 = NodeSelector(v_nodes, n2_coords)
        n2 = CircleZ(v_nodes, n2_coords, radius_disk)
        append!(nodes_BC, n2)

        #Node 3
        x_n3 = mesh_distribution[1][5]
        y_n3 = mesh_distribution[2][end-4]
        z_n3 = mesh_distribution[3][1]
        n3_coords = [x_n3, y_n3, z_n3]
        # n3 = NodeSelector(v_nodes, n3_coords)
        n3 = CircleZ(v_nodes, n3_coords, radius_disk)

        append!(nodes_BC, n3)
        #Node 4
        x_n4 = mesh_distribution[1][end-4]
        y_n4 = mesh_distribution[2][end-4]
        z_n4 = mesh_distribution[3][1]
        n4_coords = [x_n4, y_n4, z_n4]
        # n4 = NodeSelector(v_nodes, n4_coords)
        n4 = CircleZ(v_nodes, n4_coords, radius_disk)
        append!(nodes_BC, n4)
    end

    BCMat = BCMatrix(BCMat, nodes_BC, conditionBC, BC)
    # Computing the restricted degrees of freedom resulted from the applied BC
    u = zeros(num_DOF)
    (restricted_DOF, u) = RestrictedDOF(BCMat, u)

    return restricted_DOF, u
end

function _set_neumann_boundary_conditions(v_nodes, mesh_distribution, num_DOF, force_magnitude=FORCE_MAGNITUDE, outer_radius_BC=OUTER_RADIUS_BC, inner_radius_BC=INNER_RADIUS_BC, radius_F=RADIUS_BC; load_case=LOAD_CASE)
    """SETTING NODAL FORCES"""
    ForceMat = []

    x_f = mesh_distribution[1][Int((length(mesh_distribution[1]) + 1)/2)]
    y_f = mesh_distribution[2][Int((length(mesh_distribution[1]) + 1)/2)]
    z_f = mesh_distribution[3][end]

    force_coords = [x_f, y_f, z_f]   #Coordinates for the central node of the circle
    if issubset("TUBE", load_case) || load_case == "TORSION"
        force_magnitude *= -1
        # The force will be applied in a circular region
        nodes_force = CircunferenceZ(v_nodes, force_coords, inner_radius_BC, outer_radius_BC)

    elseif load_case == "SHEAR"
        nodes_force =  CircleZ(v_nodes, force_coords, 20.)

    elseif load_case == "CHAIR"
        nodes_force =  CircleZ(v_nodes, force_coords, 44.)
    end

    if load_case == "TORSION"
        # Impose the force conditions tangent to the circumference
        conditionF = [1, 1, 0] # Force applied in Z direction
        # Tangent force (resultant)
        ft = force_magnitude #/ length(nodes_force) #TODO
        # Only tangent forces → fz = 0
        fz = 0.
        # Initializes an empty array (each node is subjected to different fx, fy loads)
        f_value = [zeros(3) for _ in range(1, length=length(nodes_force), step=1)]
        # Loop in nodes_force to impose the load
        for (i, j) in enumerate(nodes_force)
            # Obtain the α angle (with the x-axis), α = atan(y_node / x_node) (sustract offset to obtain location relative to the circle)
            x_node = v_nodes[j][1] - x_f
            y_node = v_nodes[j][2] - y_f
            α = atan(y_node, x_node)
            # Assign forces
            fx = -ft * sin(α)
            fy = ft * cos(α)
            # Append to the f_value array
            f_value[i] = [fx, fy, fz]
        end

        nodes_coordinates = [vcat(v_nodes[i], j) for (i, j) in zip(nodes_force, f_value)]
        writedlm("./post-process/nodes_coordinates.csv", nodes_coordinates, ",")
    else
        conditionF = [0, 0, 1] # Force applied in Z direction

        fx1 = 0.0 / length(nodes_force)
        fy1 = 0.0 / length(nodes_force)
        fz1 = force_magnitude #/ length(nodes_force) #TODO: DELETE

        f_value = [fx1, fy1, fz1]
    end

    println("Load case: $load_case. Number of nodes: $(length(nodes_force))")
    ForceMat = BCMatrix(ForceMat, nodes_force, conditionF, f_value)
    f = zeros(num_DOF)
    (force_DOF, f) = RestrictedDOF(ForceMat, f)
    return force_DOF, f
end

# Build the mesh
_, elements, v_nodes, v_mesh, E_mesh, ν_mesh, mesh_distribution, num_DOF, nodal_connectivity, vol = build_mesh(print_paraview=true)
# Set Boundary Conditions
restricted_DOF, u, _, f = set_boundary_conditions(v_nodes, mesh_distribution, num_DOF)
# Run optimization
Enew_mesh, Enew_nodes, volfrac = run_optimization(elements, v_nodes, v_mesh, E_mesh, ν_mesh, u, f, restricted_DOF, nodal_connectivity, vol)

println("Volfrac reached f = $volfrac")

# """Pop-up window to notify that the code has been executed
# """
# win = Window("Analysis completed")
# ok = Button("Close window to finish")
# push!(win, ok)
# showall(win)

# if !isinteractive()
#     c = Condition()
#     signal_connect(win, :destroy) do widget
#         notify(c)
#     end
#     wait(c)
# end
