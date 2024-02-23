
##
"""
This paraview function is used whenever you want to visualize a mesh along with other
different results such as: strains, stresses...

ncoords must be an array containing the nodes' coordinates
mesh must be an array containing the nodes of each of the structure's elements
strains must be an array containing the data of each node related with that result (e.g. strains)
"""
function PrintParaview(ncoords, mesh, strains, stresses, displacements, E_nodes, α_nodes, α_mesh, W_nodes, name::String; verbose::Bool=true)
    #ncoords es el array que contiene las coordenadas de cada nodo
    #propiedades es el array que contiene los resultados de cada nodo: stress, strains...
    #mesh es el array que indica qué nodos pertenecen a qué elemento
    if verbose
        println("creating Paraview results file: ", name)
    end
    #println("file name ", name)

    num_pts=length(ncoords[:,1]) # 3D. #Esto es igual al numero de nodos de la estructura

    dim_tet=length(mesh[:,1]) # Tetrahedro. Esto es igual al numero de elementos de la estructura

    #a="pilar1"
    #c=".vtu"
    #nam=a+c
    nam = name

    # v_tetra should be a vector, defined previously
    #f=open("output_micro.vtu","w")
    f=open(nam,"w")
    # la funcion \n implica un salto de linea en el archivo de escritura
    write(f, """<?xml version="1.0"?>\n""")
    write(f, """<VTKFile type="UnstructuredGrid">\n""")
    write(f, "\t<UnstructuredGrid>\n")
    write(f, """\t\t<Piece NumberOfPoints=" """)
    write(f, string(Int(num_pts)))
    write(f, """ "  NumberOfCells=" """)
    write(f, string(Int(dim_tet)))
    write(f, """ ">\n""")
    write(f, "\t\t\t<Points>\n")
    write(f, """\t\t\t\t<DataArray type="Float64" NumberOfComponents="3" format="ascii">\n""")

    #Data array --> Coordinates
    #Para cada uno de los nodos (contador i) escribimos sus coordenadas
    for i in range(1,stop=num_pts)
        #coordenada x
        write(f, "\t\t\t\t\t", string(ncoords[i][1]))
        write(f, "   ")
        #coordenada y
        write(f, string(ncoords[i][2]))
        write(f, "   ")
        #coordenada z
        write(f, string(ncoords[i][3]))
        write(f, "\n")
    end

    write(f, "\t\t\t\t</DataArray>\n")
    write(f, "\t\t\t</Points>\n")
    write(f, "\t\t\t<Cells>\n")
    write(f, """\t\t\t\t<DataArray type="Int32" Name="connectivity" format="ascii">\n""")

    # Data array --> Connectivity
    # Bucle para cada uno de los elementos. Se va a seleccionar cada uno de los nodos
    # del elemento i.
    for i in range(1,stop=dim_tet)
        write(f, "\t\t\t\t\t", string(mesh[i][1]-1))
        write(f, "   ")
        write(f, string(mesh[i][2]-1))
        write(f, "   ")
        write(f, string(mesh[i][3]-1))
        write(f, "   ")
        write(f, string(mesh[i][4]-1))
        write(f, "\n")
    end

    write(f, "\t\t\t\t</DataArray>\n")
    write(f, """\t\t\t\t<DataArray type="Int32" Name="offsets" format="ascii">\n""")

    # Data array --> Offset
    #
    for i in range(1,stop=dim_tet)
        aux_p=(i)*4
        write(f, "\t\t\t\t\t", string(aux_p))
        write(f, "\n")
    end

    write(f, "\t\t\t\t</DataArray>\n")
    write(f, """\t\t\t\t<DataArray type="UInt8" Name="types" format="ascii">\n""")

    for i in range(1,stop=dim_tet)
        write(f, "\t\t\t\t\t", string(10))
        write(f, "\n")
    end

    write(f, "\t\t\t\t</DataArray>\n")
    write(f, "\t\t\t</Cells>\n")
    write(f, "\t\t\t<PointData>\n")

    # Este apartado es el que se debe copiar y pegar para cada una de las propiedades
    # que queremos estudiar.
    write(f, """\t\t\t\t<DataArray type="Float64" Name="Strains" NumberOfComponents="6" format="ascii">\n""")

    for i in range(1,stop=num_pts)
        write(f, "\t\t\t\t\t", string(strains[i][1]))
        write(f, "   ")
        write(f,  string(strains[i][2]))
        write(f, "   ")
        write(f,  string(strains[i][3]))
        write(f, "   ")
        write(f,  string(strains[i][4]))
        write(f, "   ")
        write(f,  string(strains[i][5]))
        write(f, "   ")
        write(f,  string(strains[i][6]))

        write(f, "\n")
    end

    write(f, "\t\t\t\t</DataArray>\n")
    #

    write(f, """\t\t\t\t<DataArray type="Float64" Name="Stresses" NumberOfComponents="6" format="ascii">\n""")

    for i in range(1,stop=num_pts)
        write(f, "\t\t\t\t\t", string(stresses[i][1]))
        write(f, "   ")
        write(f,  string(stresses[i][2]))
        write(f, "   ")
        write(f,  string(stresses[i][3]))
        write(f, "   ")
        write(f,  string(stresses[i][4]))
        write(f, "   ")
        write(f,  string(stresses[i][5]))
        write(f, "   ")
        write(f,  string(stresses[i][6]))

        write(f, "\n")
    end

    write(f, "\t\t\t\t</DataArray>\n")
    #

    write(f, """\t\t\t\t<DataArray type="Float64" Name="Displacements" NumberOfComponents="3" format="ascii">\n""")

    for i in range(1,stop=num_pts)
        write(f,  "\t\t\t\t\t", string(displacements[i*3-2]))
        write(f, "   ")
        write(f, string(displacements[i*3-1]))
        write(f, "   ")
        write(f, string(displacements[i*3]))

        write(f, "\n")
    end

    write(f, "\t\t\t\t</DataArray>\n")
    #

    # Young Modulus
    write(f, """\t\t\t\t<DataArray type="Float64" Name="Young Modulus" NumberOfComponents="1" format="ascii">\n""")

    for i in range(1,stop=num_pts)
        write(f, "\t\t\t\t\t", string(E_nodes[i]))

        write(f, "\n")
    end

    write(f, "\t\t\t\t</DataArray>\n")
    #

    # alphas (Young Modulus update parameter)
    write(f, """\t\t\t\t<DataArray type="Float64" Name="alpha" NumberOfComponents="1" format="ascii">\n""")

    for i in range(1,stop=num_pts)
        write(f, "\t\t\t\t\t", string(α_nodes[i]))

        write(f, "\n")
    end

    write(f, "\t\t\t\t</DataArray>\n")
    #

    #### alphas (Young Modulus update parameter) -- IN ELEMENTS INSTEAD OF NODES (NO AVERAGED)
    write(f, """\t\t\t\t<DataArray type="Float64" Name="alpha_2" NumberOfComponents="1" format="ascii">\n""")

    for i in range(1,stop=dim_tet)
        write(f, "\t\t\t\t\t", string(α_mesh[i]))

        write(f, "\n")
    end

    write(f, "\t\t\t\t</DataArray>\n")
    #

    # Elastic energy (Wi)
    write(f, """\t\t\t\t<DataArray type="Float64" Name="Elastic Energy" NumberOfComponents="1" format="ascii">\n""")

    for i in range(1,stop=num_pts)
        write(f, "\t\t\t\t\t", string(W_nodes[i]))

        write(f, "\n")
    end

    write(f, "\t\t\t\t</DataArray>\n")

    write(f, "\t\t\t</PointData>\n")
    write(f, "\t\t</Piece>\n")
    write(f, "\t</UnstructuredGrid>\n")
    write(f, "</VTKFile>\n")
    close(f)

    return;
end

##
"""
This paraview function is used whenever you want to visualize a mesh, with no stresses or
strains calculated.
ncoords must be an array containing the nodes' coordinates
mesh must be an array containing the nodes of each of the structure's elements
"""
function PrintParaview(ncoords::Array{Any,1}, mesh::Array{Any,1}, name::String; verbose::Bool=true)
    #ncoords es el array que contiene las coordenadas de cada nodo
    #propiedades es el array que contiene los resultados de cada nodo: stress, strains...
    #mesh es el array que indica qué nodos pertenecen a qué elemento
    if verbose
        println("creating Paraview mesh file: ", name)
    end

    num_pts=length(ncoords[:,1]) # 3D. #Esto es igual al numero de nodos de la estructura

    dim_tet=length(mesh[:,1]) # Tetrahedro. Esto es igual al numero de elementos de la estructura

    #a="pilar1"
    #c=".vtu"
    #nam=a*c
    nam = name

    # v_tetra should be a vector, defined previously
    #f=open("output_micro.vtu","w")
    f=open(nam,"w")
    # la funcion \n implica un salto de linea en el archivo de escritura
    write(f, """<?xml version="1.0"?> \n""")
    write(f, """<VTKFile type="UnstructuredGrid"> \n""")
    write(f, "   <UnstructuredGrid> \n")
    write(f, """      <Piece NumberOfPoints=" """)
    write(f, string(Int(num_pts)))
    write(f, """ "  NumberOfCells=" """)
    write(f, string(Int(dim_tet)))
    write(f, """ ">  \n""")
    write(f, "          <Points>  \n")
    write(f, """             <DataArray type="Float64" NumberOfComponents="3" format="ascii">   \n""")

    #Data array --> Coordinates
    #Para cada uno de los nodos (contador i) escribimos sus coordenadas
    for i in range(1,stop=num_pts)
        #coordenada x
        write(f, string(ncoords[i][1]))
        write(f, "   ")
        #coordenada y
        write(f, string(ncoords[i][2]))
        write(f, "   ")
        #coordenada z
        write(f, string(ncoords[i][3]))
        write(f, "\n")
    end

    write(f, "             </DataArray>   \n")
    write(f, "          </Points>  \n")
    write(f, "          <Cells>  \n")
    write(f, """             <DataArray type="Int32" Name="connectivity" format="ascii">   \n""")

    # Data array --> Connectivity
    # Bucle para cada uno de los elementos. Se va a seleccionar cada uno de los nodos
    # del elemento i.
    for i in range(1,stop=dim_tet)
        write(f,  string(mesh[i][1]-1))
        write(f, "   ")
        write(f, string(mesh[i][2]-1))
        write(f, "   ")
        write(f, string(mesh[i][3]-1))
        write(f, "   ")
        write(f, string(mesh[i][4]-1))
        write(f, "\n")
    end

    write(f, "              </DataArray>   \n")
    write(f, """              <DataArray type="Int32" Name="offsets" format="ascii">   \n""")

    # Data array --> Offset
    #
    for i in range(1,stop=dim_tet)
        aux_p=(i)*4
        write(f, string(aux_p))
        write(f, "\n")
    end

    write(f, "              </DataArray>   \n")
    write(f, """              <DataArray type="UInt8" Name="types" format="ascii">   \n""")

    for i in range(1,stop=dim_tet)
        write(f, string(10))
        write(f, "\n")
    end

    write(f, "              </DataArray>   \n")
    write(f, "          </Cells>  \n")

    write(f, "          <PointData>  \n")

    write(f, "          </PointData>  \n")
    write(f, "       </Piece>  \n")
    write(f, "   </UnstructuredGrid> \n")
    write(f, "</VTKFile> \n")
    close(f)

    return;
end
