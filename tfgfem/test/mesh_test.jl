using Revise
using tfgfem

sizex = 10
sizey = 1
sizez = 1
v_size = [sizex, sizey, sizez]

gapx = 1
gapy = 1
gapz = 1
v_gap = [gapx, gapy, gapz]

(nodes, elements, v_nodes, v_mesh, v_mechanical) = MeshBuild(v_size, v_gap, savedir=homedir()*"\\Documents")
