# Validation of KlocTetra and VolumeTetra functions. Results  are compared with
# example in page 246 from Book "MECANICA COMPUTACIONAL EN LA INGENIERIA CON
# APLICACIONES EN MATLAB"
using Revise
using tfgfem
#Defining the tethaedral element

v_nodes = []
push!(v_nodes, [0.0, 0.0, 0.0])
push!(v_nodes, [1.0, 0.5, 0.5])
push!(v_nodes, [0.0, 0.5, 1.0])
push!(v_nodes, [0.0, 0.0, 1.0])

v_mesh = []
push!(v_mesh, [1,2,3,4])

nodes = length(v_nodes)
elements = length(v_mesh)

# Obtaining the nodal connectivity
nodal_connectivity = NodalConnectivity(v_nodes, v_mesh)

# Obtaining the elemental connectivity
(coord_glob_element, connect_element_glob) = ElementConnectivity(v_nodes, v_mesh)

E = 20000
pois = 0.2

# Initializing mechanical properties in an array
E_mesh = []
push!(E_mesh, E)
pois_mesh = []
push!(pois_mesh, pois)

gdl_node = 3
ngdl = nodes*gdl_node

aux_vec_coord = []
for j2 in range(1,stop=4)
    #En naux se almacena
    naux=v_mesh[1][j2]
    for k2 in range(1,stop=3)  #recorre las 3 columnas de v_nodes
        push!(aux_vec_coord, v_nodes[naux][k2])
    end
end

voltet = VolumeTetra(aux_vec_coord)

(Kloc, B, voltet, D) = KlocTetra(aux_vec_coord, E, pois)
Kglob = Kloc * voltet

#vector fuerzas nodales
f = [0,0,0,1851.85,0,0,0,0,0,0,0,0]

# nodes 1,3 and 4 are embedded (dsiplacements = 0)
# u = [0,0,0,u4,u5,u6,0,0,0,0,0,0]
# u initialized as zeros
u = zeros(ngdl)
free_DOF = [4,5,6]
restrict_DOF = setdiff(1:ngdl, free_DOF)

(Kglob_FF, Kglob_FR, Kglob_RR, u_F, u_R, u, f_F, f_R) = SolveFEM(Kglob, u, f, restrict_DOF)

u_F = Kglob_FF\f_F

for idof = 1:length(free_DOF)
    DOF = free_DOF[idof]
    u[DOF] = u_F[idof]
end

for idof = 1:length(restrict_DOF)
    DOF = restrict_DOF[idof]
    u[DOF] = u_R[idof]
end


(strains_el, stresses_el) = StrainsStressesElem(v_nodes, v_mesh, E_mesh, pois_mesh, u)
(voigtstrains_el, voigtstresses_el) = VoigtFormula(strains_el, stresses_el)

(strains_nod, stresses_nod) = StrainsStressesNodes(nodal_connectivity, strains_el, stresses_el)
