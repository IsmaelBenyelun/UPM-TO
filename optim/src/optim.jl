module optim

# Loading dependencies
using LinearAlgebra
using DelimitedFiles
using tfgfem

include("DrMatrix.jl")
include("invDrMatrix.jl")
export DrMatrix
export invDrMatrix

include("CubeLayering.jl")
export ElementLayer
export CubeLayeringZ

include("ConstantH.jl")
export ConstantH

end # module
