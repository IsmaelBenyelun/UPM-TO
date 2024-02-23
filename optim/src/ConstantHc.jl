function ConstantHc(invDr::Array{Float64,2}, voigtstresses::Array{Float64,1}, vol_el::Float64)
"""ComplementaryGD add-on W^c_i = H^c / E_i
"""
    if length(voigtstresses)==6
        aux = invDr*voigtstresses
        aux2 = transpose(aux)*voigtstresses
        Hc = aux2*vol_el*0.5

    elseif length(voigtstresses)==3
        invDr = invDr[1:3,1:3]
        aux = invDr*voigtstresses
        aux2 = transpose(aux)*voigtstresses
        Hc = aux2*vol_el*0.5

    elseif length(voigtstresses)==1
        invDr = invDr[1:1,1:1]
        aux = invDr*voigtstresses
        aux2 = transpose(aux)*voigtstresses
        Hc = aux2*vol_el*0.5
    end
    return Hc
end
