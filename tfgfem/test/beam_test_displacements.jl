using Revise
using tfgfem

cd("C:\\Users\\Hugo\\Desktop\\tfgfem\\Beam Validation")

"""MESH BUILDING"""
# Meshing declaration
sizex = 1000
sizey = 50
sizez = 50
v_size = [sizex, sizey, sizez]

gapx = 10.0
gapy = 10.0
gapz = 10.0
v_gap = [gapx, gapy, gapz]

jumpx = 0

# Mechanical properties
E=1000.0
pois=0.3

# Obtaining vector por nodes, mesh and mechanical properties
(nodes, elements, v_nodes, v_mesh, v_mechanical) = MeshBuild(v_size, v_gap, E, pois)
(E_mesh, pois_mesh) = MechanicalMesh(v_mechanical)
mesh_distribution = MeshDistribution(v_size, v_gap, jumpx)
# Printing the obtained mesh in paraview
PrintParaview(v_nodes, v_mesh, "beam_mesh.vtu")

dof_node = 3
num_DOF = nodes*dof_node


"""CONNECTIVITY"""
# Obtaining the nodal connectivity
#nodal_connectivity = NodalConnectivity(v_nodes, v_mesh)

# Obtaining the elemental connectivity
#(element_info, connect_element_glob) = ElementConnectivity(v_nodes, v_mesh)


"""BOUNDARY CONDITIONS"""
#Position and condition
dimension1 = [1,0,0]
pos_BC1 = [-240.0, 0.0, 0.0]
restr_displacements1 = [1,1,1]
BC1 = [0.0, 0.0, 0.0]
nodes_BC1 = NodeGroup(v_nodes, dimension1, pos_BC1)

dimension2 = [1,0,1]
pos_BC2 = [260.0, 0.0, 60.0]
restr_displacements2 = [0,0,1]
BC2 = [0.0, 0.0, -5.0]
nodes_BC2 = NodeGroup(v_nodes, dimension2, pos_BC2)

# Nodes where each BC is applied, restricted degrees of freedom and displacements (u)
#(nodes_BC11, nodes_BC22, restrict_DOF, u) = BoundaryConditions(v_nodes, elements, pos_BC1, pos_BC2, BC1, BC2)

#Computing Boundary Condition Matriz, BCMat.
BCMat = []
BCMat = BCMatrix(BCMat, nodes_BC1, restr_displacements1, BC1)
BCMat = BCMatrix(BCMat, nodes_BC2, restr_displacements2, BC2)

u = zeros(num_DOF)
(restricted_DOF, u) = RestrictedDOF(BCMat, u)

f = zeros(num_DOF)

"""RUNNING THE FINITE ELEMENT ANALYSIS"""
(voigtstrains_el, voigtstresses_el, strains_nod, stresses_nod, strains_el, stresses_el) =
FEAnalysis(v_nodes, v_mesh, E_mesh, pois_mesh, u, f, restricted_DOF)
nodal_connectivity = NodalConnectivity(v_nodes, v_mesh)
E_nodes = NodalYoungModule(nodal_connectivity, E_mesh)
PrintParaview(v_nodes, v_mesh, strains_nod, stresses_nod, u, E_nodes, "beam_analysis_displacements.vtu")
