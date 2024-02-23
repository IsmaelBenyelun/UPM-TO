using Revise
using tfgfem

cd("C:\\Users\\Hugo\\Desktop\\tfgfem")

# Meshing declaration
sizex = 10
sizey = 10
sizez = 10
v_size = [sizex, sizey, sizez]

gapx = 1.0
gapy = 1.0
gapz = 1.0
v_gap = [gapx, gapy, gapz]

jumpx = 0

# Mechanical properties
E=28300.0
pois=0.2

# Obtaining vector por nodes, mesh and mechanical properties
(nodes, elements, v_nodes, v_mesh, v_mechanical) = MeshBuild(v_size, v_gap, E, pois)
mesh_distribution = MeshDistribution(v_size, v_gap, jumpx)


# Printing the obtained mesh in paraview
PrintParaview(v_nodes, v_mesh)

""" CONNECTIVITY """

# Obtaining the nodal connectivity
nodal_connectivity = NodalConnectivity(v_nodes, v_mesh)

# Obtaining the elemental connectivity
(element_info, connect_element_glob) = ElementConnectivity(v_nodes, v_mesh)


# Initializing mechanical properties in an array
(E_mesh, pois_mesh) = MechanicalMesh(v_mechanical)

""" COMPUTING GLOBAL K MATRIX """
gdl_node = 3
(Kglob, collect_aux_vec_coord, B) = Kglobal(gdl_node, v_nodes, v_mesh, E_mesh, pois_mesh)
num_DOF = length(Kglob[:,1])
##
""""SETTING BOUNDARY CONDITIONS"""

# BOUNDARY CONDITIONS 2. EMBEDDED
condition2 = [1,1,1]
BC2 = [0.0,0.0,0.0]

# Selecting affected nodes
x_s2 = mesh_distribution[1][Int((length(mesh_distribution[1]) + 1)/2)]
y_s2 = mesh_distribution[2][Int((length(mesh_distribution[1]) + 1)/2)]
z_s2 = mesh_distribution[3][1]
node_coords2 = [x_s2, y_s2, z_s2]
center_node2 = NodeSelector(v_nodes, node_coords2)
radius2 = 1.0
# CircleZ will be the vector containing the nodes affected by the future BC.
nodes_BC2 = CircleZ(v_nodes, node_coords2, radius2)

BCMat = []
BCMat = BCMatrix(BCMat, nodes_BC2, condition2, BC2)

# Computing the restricted degrees of freedom resulted from the applied BC
u = zeros(num_DOF)
(restricted_DOF, u) = RestrictedDOF(BCMat, u)


##
"""SETTING NODAL FORCES"""
condition1 = [0,0,1]
fx1 = 0.0
fy1 = 0.0
fz1 = 1350.0
f_value1 = [fx1, fy1, fz1]

x_s1 = mesh_distribution[1][Int((length(mesh_distribution[1]) + 1)/2)]
y_s1 = mesh_distribution[2][Int((length(mesh_distribution[1]) + 1)/2)]
z_s1 = mesh_distribution[3][end]
node_coords1 = [x_s1, y_s1, z_s1]

center_node1 = NodeSelector(v_nodes, node_coords1)
radius1 = 2.0

nodes_force1 =  CircleZ(v_nodes, node_coords1, radius1)

ForceMat = []
ForceMat = BCMatrix(ForceMat, nodes_force1, condition1, f_value1)

f = zeros(num_DOF)
(force_DOF, f) = RestrictedDOF(ForceMat, f)

##
"""SOLVING THE FEM MODEL"""

(Kglob_FF, Kglob_FR, Kglob_RR, u_F, u_R, u, f_F, f_R) = SolveFEM(Kglob, u, f, restricted_DOF)


##
"""COMPUTING STRAINS AND STRESSES"""
# Strains and stresses inside the element
(strains_el, stresses_el) = StrainsStressesElem(v_nodes, v_mesh, E_mesh, pois_mesh, u)
(voigtstrains_el, voigtstresses_el) = VoigtFormula(strains_el, stresses_el)
# Strains and stresses calculated in each node
(strains_nod, stresses_nod) = StrainsStressesNodes(nodal_connectivity, strains_el, stresses_el)

E_nodes = NodalYoungModule(nodal_connectivity, E_mesh)
PrintParaview(v_nodes, strains_nod, stresses_nod, u, E_nodes, v_mesh)
