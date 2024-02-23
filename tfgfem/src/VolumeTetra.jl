# Volume calculation for tethraedral elements.
function VolumeTetra(xev,yev,zev,xn1,yn1,zn1,xn2,yn2,zn2,xn3,yn3,zn3)
    #=global mat

    x1=1.0*(xn1-xev)
    y1=1.0*(yn1-yev)
    z1=1.0*(zn1-zev)

    x2=1.0*(xn2-xev)
    y2=1.0*(yn2-yev)
    z2=1.0*(zn2-zev)

    x3=1.0*(xn3-xev)
    y3=1.0*(yn3-yev)
    z3=1.0*(zn3-zev)

    #mat=[[x1 x2 x3];[y1 y2 y3];[z1 z2 z3]]
    mat=[[x1 y1 z1];[x2 y2 z2];[x3 y3 z3]]
    =#

    (x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4) = (xev,yev,zev,xn1,yn1,zn1,xn2,yn2,zn2,xn3,yn3,zn3)
    mat = [[x1 y1 z1 1];[x2 y2 z2 1];[x3 y3 z3 1];[x4 y4 z4 1]]


#=
    if det(mat)<0
        print("negative determinant ")
    elseif det(mat)==0
        print("null determinant ")
    else
        print("positive determinant ")
    end
=#

    volume= (1.0/6.0)*(det(mat))


    return (volume);

end

##
function VolumeTetra(coords_element::Array{Any,1})

    #(xev,yev,zev,xn1,yn1,zn1,xn2,yn2,zn2,xn3,yn3,zn3) = coords_element

    #=
    x1=1.0*(xn1-xev)
    y1=1.0*(yn1-yev)
    z1=1.0*(zn1-zev)

    x2=1.0*(xn2-xev)
    y2=1.0*(yn2-yev)
    z2=1.0*(zn2-zev)

    x3=1.0*(xn3-xev)
    y3=1.0*(yn3-yev)
    z3=1.0*(zn3-zev)

    #mat=[[x1 x2 x3];[y1 y2 y3];[z1 z2 z3]]
    mat=[[x1 y1 z1];[x2 y2 z2];[x3 y3 z3]]
    =#

    (x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4) = coords_element

    mat = [[x1 y1 z1 1];[x2 y2 z2 1];[x3 y3 z3 1];[x4 y4 z4 1]]

#=
    if det(mat)<0
        print("negative determinant")
    elseif det(mat)==0
        print("null determinant")
    else
        print("positive determinant")
    end
=#

    volume= (1.0/6.0)*(det(mat))


    return (volume);

end
