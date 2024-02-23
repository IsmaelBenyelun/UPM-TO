using Revise
using tfgfem
using optim
using TickTock

cd("C:\\Users\\Hugo\\Desktop\\optim\\LOAD CASE 1")
LC=1
println("Running cube analysis for further optimization")

"""MESH BUILDING"""
# Meshing declaration.
sizex = 100
sizey = 100
sizez = 100
v_size = [sizex, sizey, sizez]

gapx = sizex/30.0
gapy = sizex/30.0
gapz = sizex/30.0
v_gap = [gapx, gapy, gapz]

jumpx = 0

# Mechanical properties
E=28300.0
pois=0.2

# Obtaining vector por nodes, mesh and mechanical properties
(nodes, elements, v_nodes, v_mesh, v_mechanical) = MeshBuild(v_size, v_gap, E, pois)
(E_mesh, pois_mesh) = MechanicalMesh(v_mechanical)
mesh_distribution = MeshDistribution(v_size, v_gap, jumpx)
PrintParaview(v_nodes, v_mesh, "mesh.vtu")

dof_node = 3
num_DOF = nodes*dof_node

@time (Kglob, B) = Kglobal(dof_node, v_nodes, v_mesh, E_mesh, pois_mesh)
@time (KglobSp, B) = KglobalSp(dof_node, v_nodes, v_mesh, E_mesh, pois_mesh)

##
"""SETTING BOUNDARY CONDITIONS"""

# BOUNDARY CONDITIONS. EMBEDDED CIRCLE AT THE BOTTOM OF THE CUBE
conditionBC = [1,1,1]    #All the displacements are restricted
BC = [0.0,0.0,0.0]     #Value of the restrictions

# Selecting affected nodes
x_bc = mesh_distribution[1][Int((length(mesh_distribution[1]) + 1)/2)]
y_bc = mesh_distribution[2][Int((length(mesh_distribution[1]) + 1)/2)]
z_bc = mesh_distribution[3][1]
bc_coords = [x_bc, y_bc, z_bc]

radiusBC = 20.0
# CircleZ will be the vector containing the nodes affected by the future BC.
nodes_BC = CircleZ(v_nodes, bc_coords, radiusBC)

BCMat = []
BCMat = BCMatrix(BCMat, nodes_BC, conditionBC, BC)

# Computing the restricted degrees of freedom resulted from the applied BC
u = zeros(num_DOF)
(restricted_DOF, u) = RestrictedDOF(BCMat, u)


##
"""SETTING NODAL FORCES"""
# The force will be applied in a circular region
conditionF = [0,0,1]   #Force applied in Z direction
fx1 = 0.0
fy1 = 0.0
fz1 = 1000.0
f_value = [fx1, fy1, fz1]

x_f = mesh_distribution[1][Int((length(mesh_distribution[1]) + 1)/2)]
y_f = mesh_distribution[2][Int((length(mesh_distribution[1]) + 1)/2)]
z_f = mesh_distribution[3][end]
force_coords = [x_f, y_f, z_f]   #Coordinates for the central node of the circle

radiusF = 20.0

nodes_force =  CircleZ(v_nodes, force_coords, radiusF)

ForceMat = []
ForceMat = BCMatrix(ForceMat, nodes_force, conditionF, f_value)

f = zeros(num_DOF)
(force_DOF, f) = RestrictedDOF(ForceMat, f)

"""RUNNING THE FINITE ELEMENT ANALYSIS"""

@time (voigtstrains_el, voigtstresses_el, strains_nod, stresses_nod, strains_el, stresses_el) = FEAnalysis(v_nodes,
v_mesh, E_mesh, pois_mesh, u, f, restricted_DOF)

nodal_connectivity = NodalConnectivity(v_nodes, v_mesh)
E_nodes = NodalYoungModule(nodal_connectivity, E_mesh)
PrintParaview(v_nodes, v_mesh, strains_nod, stresses_nod, u, E_nodes, "Cube_analysis_LC1.vtu")

(K, B) = Kglobal(dof_node, v_nodes, v_mesh, E_mesh, pois_mesh)
uT = transpose(u)
compliance = 0.5 * uT * K * u
