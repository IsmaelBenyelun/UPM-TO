using SparseArrays

"""
    (Kloc, Bmat, voltet, D, straMat) = KlocTetra(aux_vec_coord, E_el, ν_el)

Compute the element K matrix (Kloc) given a coordinates vector and the mechanical
properties
Bmat: matrix including shape functions' derivative terms
voltet: element volume
D: constitutive matrix
straMat: D*Bmat
"""
# calculation of the local K matrix in tethaedral elements
function KlocTetra(aux_vec_coord, E_el, ν_el)

    voltet = VolumeTetra(aux_vec_coord)

    Bmat = BMatrix(aux_vec_coord)
    #Bmat=sparse(Bmat_aux) #Sparse transformations to obtain a reducction in computing time, is a vector of vectors so it´s bound not to work

    Bt = transpose(Bmat)
    #Bt=sparse(Bt_aux) This is the try of computing B with sparse

    D_aux = DMatrix(E_el, ν_el)
    D=sparse(D_aux)


    # Compute K
    M1=Bt*D
    Kloc_aux=M1*Bmat
    Kloc=sparse(Kloc_aux)
    straMat=D*Bmat

    return (Kloc, Bmat, voltet, D, straMat)

end


"""
    B = BMatrix(Bloc)

Compute the matrix including the derivative of the shape figures.
The matrix is computed in order to match the strains formulation:
(εx, εy, εz, γxy, γyz, γxz)

aux_vec_coord: array including the 12 coordinates (3 per node) of a tethaedron

"""


function BMatrix(aux_vec_coord)

    F = zeros(4,4)
    P = zeros(4,4)

    ## Definition of transposition of Jacobien Matrix
    #=La matriz F es igual a la matriz V definida en el libro de FEM=#
    for it in range(1,stop=4)
        F[it,:] = [1 aux_vec_coord[(it-1)*3+1] aux_vec_coord[(it-1)*3+2] aux_vec_coord[(it-1)*3+3]]
        #println("row", F[it,:])
    end

    if det(F)==0
        println("F error", F, det(F))
        return (0)
    end

    Bloc=[]
    ## Calculation of each Derivative Terms
    for it in range(1,stop=3)
        push!(Bloc, [])
        for j in range(1,stop=4)
            Pmat = copy(F)
            #=eliminamos columnas. Esta operacion determina si calculamos el termino
             ai, bi, ci o hi. Como en la formula de la matriz Bi (página 244 libros Chaves)
             no aparece el término ai, el contador de columnas empieza por it+1=#
            Pmat = Pmat[:, 1:end .!=(it+1)]
            #=eliminamos filas. Esta operacion determina el subindices de las letras.
            Es decir, a1, b1,c1 y h1 (j=1); a2, b2, c2, h2 (j=2):=#
            Pmat = Pmat[1:end .!=j, :]
            if (it == j)  ||  ((it == 3) && (j == 1))  ||  ((it == 1) && (j == 3))  ||  ((it == 2) && (j == 4))
                push!(Bloc[it], -det(Pmat)/det(F))
            else
                push!(Bloc[it], det(Pmat)/det(F))
            end
        end
    end

    # B declaration
    B=zeros(6, 3*4)
    #Estas 3 primeras filas corresponden  al calculo de εx εy εz.
    B[1,:]=[Bloc[1][1] 0 0 Bloc[1][2] 0 0 Bloc[1][3] 0 0 Bloc[1][4] 0 0]
    B[2,:]=[0 Bloc[2][1] 0 0 Bloc[2][2] 0 0 Bloc[2][3] 0 0 Bloc[2][4] 0]
    B[3,:]=[0 0 Bloc[3][1] 0 0 Bloc[3][2] 0 0 Bloc[3][3] 0 0 Bloc[3][4]]

    #Estas 3 filas corresponden con el cálculo de 0.5γxy 0.5γyz 0.5γxz respectivamente
    # B[4,:]=[0.5*Bloc[2][1] 0.5*Bloc[1][1] 0 0.5*Bloc[2][2] 0.5*Bloc[1][2] 0 0.5*Bloc[2][3] 0.5*Bloc[1][3] 0 0.5*Bloc[2][4] 0.5*Bloc[1][4] 0]
    # B[6,:]=[0.5*Bloc[3][1] 0 0.5*Bloc[1][1] 0.5*Bloc[3][2] 0 0.5*Bloc[1][2] 0.5*Bloc[3][3] 0 0.5*Bloc[1][3] 0.5*Bloc[3][4] 0 0.5*Bloc[1][4]]
    # B[5,:]=[0 0.5*Bloc[3][1] 0.5*Bloc[2][1] 0 0.5*Bloc[3][2] 0.5*Bloc[2][2] 0 0.5*Bloc[3][3] 0.5*Bloc[2][3] 0 0.5*Bloc[3][4] 0.5*Bloc[2][4]]

    # Second option: without multiplying by 0.5 last 3 rows. This means: γxy γyz γxz
    B[4,:]=[Bloc[2][1] Bloc[1][1] 0 Bloc[2][2] Bloc[1][2] 0 Bloc[2][3] Bloc[1][3] 0 Bloc[2][4] Bloc[1][4] 0]
    B[5,:]=[0 Bloc[3][1] Bloc[2][1] 0 Bloc[3][2] Bloc[2][2] 0 Bloc[3][3] Bloc[2][3] 0 Bloc[3][4] Bloc[2][4]]
    B[6,:]=[Bloc[3][1] 0 Bloc[1][1] Bloc[3][2] 0 Bloc[1][2] Bloc[3][3] 0 Bloc[1][3] Bloc[3][4] 0 Bloc[1][4]]



    return B

end


"""
    D = DMatrix(E_el, ν_el)

Compute the constitutive matrix D
"""

function DMatrix(E_el, ν_el)

    #Calculamos las constantes de lamé: λ y μ. La fórmula está en la pag 245 de Chaves
    λ = ν_el*E_el / ((1 + ν_el) * (1 - 2*ν_el))
    μ = E_el/(2*(1 + ν_el))

    #Estos son los diferentes terminos que forman la matriz constitutiva clásica, D
    d1 = λ + 2*μ
    d2 = λ
    d3 = μ
    # d3 = 2.0*μ # Wrong!!

    # Compute D
    D = [d1 d2 d2 0 0 0; d2 d1 d2 0 0 0; d2 d2 d1 0 0 0; 0 0 0 d3 0 0; 0 0 0 0 d3 0; 0 0 0 0 0 d3]

    return D

end
