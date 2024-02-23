module tfgfem

# Loading dependencies
using LinearAlgebra
using DelimitedFiles
using IterativeSolvers
using SparseArrays

    # including function for printing in paraview
    include("PrintParaview.jl")
    export PrintParaview

    # including function for mesh building
    include("MeshBuild.jl")
    export MeshBuild
    export MechanicalMesh
    export NodalYoungModule

    include("MeshDistribution.jl")
    export MeshDistribution

    include("MeshFilter.jl")
    export MeshFilter

    include("NodeSelection.jl")
    export CircleZ
    export Circle
    export CircunferenceZ
    export NodeSelector
    export NodeGroup
    export NodeLayerZ

    # function for node connectivity
    include("Connectivity.jl")
    export ElementConnectivity
    export NodalConnectivity
    export CoordsElement
    export ElementCenterCoord

    # function for calculating volumnes in tethaedros
    include("VolumeTetra.jl")
    export VolumeTetra

    # calculation of the local K matrix in tethraedral elements
    include("KlocTetra.jl")
    export KlocTetra
    export BMatrix
    export DMatrix

    # calculation of the global K matrix of the structure
    include("Kglobal.jl")
    export Kglobal
    export KglobalSp

    # turn BC into vectors for further application
    include("BoundaryConditions.jl")
    export BoundaryConditions
    export BCMatrix

    include("RestrictedDOF.jl")
    export RestrictedDOF


    include("KglobalLL.jl")
    export KglobalLL

    # BC application and solving
    include("SolveFEM.jl")
    export SolveFEM

    # calculation of the reactions

    # Strains, stresses... calculation in elements
    include("StrainsStresses.jl")
    export StrainsStressesElem
    export StrainsStressesNodes
    export VoigtFormula


    include("FEAnalysis.jl")
    export FEAnalysis

end  # module
