""" CASE 2 PLOTS """
#aux2 VALUE PLOTTING
for n = 1:10
    display(scatter(aux_collector[n], title = "aux values for loop $n", ylims = (-1,10)))
    #savefig("C:\\Users\\Hugo\\Desktop\\optim\\aux2 plots\\aux2 values$n.png")
end

# Hi VALUE PLOTTING
for n = 1:10
    scatter(Hi_collector[n], title = "Hi values for loop $n")#, ylims = (-0.0000031,0.0000031))
    meanH = mean(Hi_collector[n])
    display(hline!([meanH]))
    #savefig("C:\\Users\\Hugo\\Desktop\\optim\\Hi2 plots\\Hi2 values$n.png")
end

# Ei2 VALUE PLOTTING
for n = 1:10
    display(scatter(Ei2collector[n], title = "Ei2 values for loop $n"))#, ylims = (-1e5,1e7)))
    #savefig("C:\\Users\\Hugo\\Desktop\\optim\\Layer 10 Analysis\\Ei2 plots\\Stiffness problem\\Ei2 zoom values$n.png")
end

#Wi2 VALUE PLOTTING
for n = 1:15
    display(scatter(Wi2collector[n], title = "Wi2 values for loop $n", ylims = (0,10.0)))
    #savefig("C:\\Users\\Hugo\\Desktop\\optim\\Wi2 plots\\Wi2 values$n.png")
end

#COMPLIANCE
plot(workcollector2, title = "layer work", ylims = (0,400))
savefig("C:\\Users\\Hugo\\Desktop\\optim\\compliance when switching negative values for their absolute value.png")

plot(compliance_collector, title = "Compliance evolution. LC3 (V5.3) (Emax=$Emax)", ylims = (0,3000))
savefig("C:\\Users\\Hugo\\Desktop\\optim\\Compliance evolution. LC3 (V5.2) (Emax=$Emax).png")


""" COMPARISON """
#aux COMPARISON
for n = 1:10
    scatter(auxcollector[n], title = "aux values comparison loop $n",
    label = "aux", ylims = (-1.0,10.0))
    display(scatter!(aux2collector[n], label = "aux2"))
    savefig("C:\\Users\\Hugo\\Desktop\\optim\\comparison aux1 and aux2\\aux plots\\comp_aux_$n.png")
end


# Hi COMPARISON
for n = 1:10
    scatter(Hicollector[n], title = "Hi values comparison loop $n",
    label = "Hi", ylims = (-0.00031,0.00031))
    meanH = mean(Hicollector[n])
    hline!([meanH], label = "meanHi")
    scatter!(Hi2collector[n], label = "Hi2")
    meanH2 = mean(Hi2collector[n])
    display(hline!([meanH2], label = "meanHi2"))
    savefig("C:\\Users\\Hugo\\Desktop\\optim\\comparison aux1 and aux2\\Hi plots\\comp_Hi_$n.png")

end

# Ei COMPARISON
for n = 1:10
    scatter(Eicollector[n], title = "Ei values for loop $n", label = "Ei", ylims = (-1e7,3e7))
    display(scatter!(Ei2collector[n], label = "Ei2"))
    savefig("C:\\Users\\Hugo\\Desktop\\optim\\comparison aux1 and aux2\\Ei plots\\comp_Ei_$n.png")
end

#Wi COMPARISON
for n = 1:10
    scatter(Wicollector[n], label = "Wi", title = "Wi values for loop $n", ylims = (0,9.0))
    display(scatter!(Wi2collector[n], label = "Wi2"))
    savefig("C:\\Users\\Hugo\\Desktop\\optim\\comparison aux1 and aux2\\Wi plots\\comp_Wi_$n.png")
end

#WORK COMPARISON

plotcomparison = scatter(workcollector, label = "compliance", title = "compliance evolution comparison")
scatter!(plotcomparison, workcollector2, label = "compliance2")
savefig("C:\\Users\\Hugo\\Desktop\\optim\\comparison aux1 and aux2\\compliance_comparison30.png")


plot(auxcollector, title="Evolution aux value for each element layer")
