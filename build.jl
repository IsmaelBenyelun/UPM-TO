import Pkg
tfgfem_path = joinpath(pwd(), "tfgfem")
optim_path = joinpath(pwd(), "optim")
Pkg.develop(path=tfgfem_path)
Pkg.develop(path=optim_path)

print("Testing packages... ")
using tfgfem
using optim
println("Ok!")

print("Testing module... ")
push!(LOAD_PATH, pwd())
using Optimfunction
println("Ok!")
