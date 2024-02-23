using Revise
using tfgfem

cd("C:\\Users\\Hugo\\Desktop\\tfgfem\\Beam Validation")

"""MESH BUILDING"""
# Meshing declaration
sizex = 500
sizey = 50
sizez = 50
v_size = [sizex, sizey, sizez]

gapx = sizex/100.0
gapy = sizey/10.0
gapz = sizez/10.0
v_gap = [gapx, gapy, gapz]

jumpx = 0

# Mechanical properties
E=30000.0
I = sizey^4/12
pois=0.3
G = E/(2*(1+pois))
A = sizez * sizey


# Obtaining vector por nodes, mesh and mechanical properties
(nodes, elements, v_nodes, v_mesh, v_mechanical) = MeshBuild(v_size, v_gap, E, pois)
(E_mesh, pois_mesh) = MechanicalMesh(v_mechanical)
mesh_distribution = MeshDistribution(v_size, v_gap, jumpx)
# Printing the obtained mesh in paraview
PrintParaview(v_nodes, v_mesh, "beam_mesh.vtu")

dof_node = 3
num_DOF = nodes*dof_node


"""BOUNDARY CONDITIONS"""
#Position and condition
dimension1 = [1,0,0]
pos_BC1 = [mesh_distribution[1][1], 0.0, 0.0]
restr_displacements1 = [1,1,1]
BC1 = [0.0, 0.0, 0.0]
nodes_BC1 = NodeGroup(v_nodes, dimension1, pos_BC1)

#Computing Boundary Condition Matriz, BCMat.
BCMat = []
BCMat = BCMatrix(BCMat, nodes_BC1, restr_displacements1, BC1)

u = zeros(num_DOF)
(restricted_DOF, u) = RestrictedDOF(BCMat, u)


"""SETTING NODAL FORCES"""
# The force will be applied in a circular region
conditionF = [0,0,1]   #Force applied in Z direction
fx1 = 0.0
fy1 = 0.0
fz1 = -1000.0
f_value = [fx1, fy1, fz1]

dimension2 = [1,0,1]
pos_forces = [mesh_distribution[1][end], 0.0, mesh_distribution[3][end]]
nodes_force = NodeGroup(v_nodes, dimension2, pos_forces)

ForceMat = []
ForceMat = BCMatrix(ForceMat, nodes_force, conditionF, f_value)

f = zeros(num_DOF)
(force_DOF, f) = RestrictedDOF(ForceMat, f)


"""RUNNING THE FINITE ELEMENT ANALYSIS"""
(voigtstrains_el, voigtstresses_el, strains_nod, stresses_nod, strains_el, stresses_el) =
FEAnalysis(v_nodes, v_mesh, E_mesh, pois_mesh, u, f, restricted_DOF)
nodal_connectivity = NodalConnectivity(v_nodes, v_mesh)
E_nodes = NodalYoungModule(nodal_connectivity, E_mesh)
PrintParaview(v_nodes, v_mesh, strains_nod, stresses_nod, u, E_nodes, "beam_analysis_forces.vtu")

max_displ = maximum(u)
min_displ = minimum(u)

w_Bern = (fz1*sizex^3)/(3*E*I)

w_Timoshenko = w_Bern + (fz1*sizex)/(G*5*A/6)

L = sizex
H = sizey
x = L
y = H
P = fz1
w_Timo_Good = P/(6*E*I) * (3*pois*(y-H/2)^2*(L-x) + 1/4*(4+5*pois)*H^2*x + (3*L-x)*x^2)
