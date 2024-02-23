

function ConstantH(Dr::Array{Float64,2}, voigtstrains::Array{Float64,1}, vol_el::Float64)


    if length(voigtstrains)==6
        aux = Dr*voigtstrains
        aux2 = transpose(aux)*voigtstrains
        H = aux2*vol_el*0.5

    elseif length(voigtstrains)==3
        Dr = Dr[1:3,1:3]
        aux = Dr*voigtstrains
        aux2 = transpose(aux)*voigtstrains
        H = aux2*vol_el*0.5

    elseif length(voigtstrains)==1
        Dr = Dr[1:1,1:1]
        aux = Dr*voigtstrains
        aux2 = transpose(aux)*voigtstrains
        H = aux2*vol_el*0.5
    end
    return H
end
