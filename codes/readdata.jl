using HDF5
using Plots
using LaTeXStrings
using Measures
#=
plot(kpoints[3,:], Eks,
     #markershape=:circle,
     xlabel=L"k_z/2\pi",
     ylabel=L"E(\mathbf{k})",
	framestyle = :box,grid=false,
    xtickfontsize=12, ytickfontsize=12,
    xguidefontsize=12, yguidefontsize=12,
    legendfontsize=10,
    legend = false,
    margin=2mm,
    )
#hline!([-1.143,-1.124,-0.8,-0.352,-0.312,0,0.312,0.352,0.8,1.124,1.143],linestyle=:dash)
hline!([-0.4,0.4],color=:black,linestyle=:dash)
#vline!([-0.318,-0.187,0.187,0.318], linestyle=:dash)
#vline!([-0.25,0.25],color=:black, linestyle=:dash)
savefig("~/Desktop/dispersionGam08.pdf")
=#
#=
f = h5open("EDcpgeL16OBCO0sigman2.h5","r")
xs = read(f, "xs")
cpgeL16 = read(f, "cpge")
w1 = read(f, "w1")
w2 = read(f, "w2")
close(f)

f = h5open("EDcpgeL8TBC100O0sigman2.h5","r")
xs = read(f, "xs")
cpges = read(f, "cpges")
w1 = read(f, "w1")
w2 = read(f, "w2")
close(f)
cpgeML8 = sum(cpges,dims=2)./100

f = h5open("EDcpgeL14OBCO0sigman2.h5","r")
cpgeL14 = read(f, "cpge")
close(f)

f = h5open("EDcpgeL8OBCO0sigman2.h5","r")
cpgeL8 = read(f, "cpge")
close(f)

f = h5open("EDcpgeL10OBCO0sigman2.h5","r")
cpgeL10 = read(f, "cpge")
close(f)

@show factor = maximum(cpgeL16)./maximum(cpgeML8)
#@show factor = sum(abs.(cpgeL16))./sum(abs.(cpgeML8))
plot(xs, cpgeML8,
     label = L"L=8, TBC 100",
     xlabel = L"E_F/D",
     ylabel = "CPGE",
     xlims = [-0.4,0.4],
     #title = L"L=8, TBC 100",
    margin=3mm,)
plot!(xs, cpgeL8,
     label = L"L=8, PBC")
plot!(xs, cpgeL10,
     label = L"L=10, PBC")
#plot!(xs, cpgeL14,label = L"L=14, PBC")
plot!(xs, cpgeL16,
     label = L"L=16, PBC")
vline!([-1.6,-0.8,0,0.8,1.6]./5.2,linestyle=:dash,
      label = false)
=#
#=
f = h5open("KPMcpgeL100NCs.h5","r")
Efs = read(f, "Efs")
cpgeKPMNCs = read(f, "CPGEs")
close(f)

plot(xs, cpgeML8./600,
     label = L"L=8, TBC 100, max=600",
     xlabel = L"E_F/D",
     ylabel = "CPGE",
     xlims = [-0.4,0.4],
     #title = L"L=8, TBC 100",
    margin=3mm,)
plot!(xs, cpgeL16./8150,
     label = L"L=16, OBC, max=8150")
plot!(Efs, imag.(cpgeKPMNCs[:,5])./1403,
     label = L"L=100, N_C=256, max=1403")
vline!([-1.6,-0.8,0,0.8,1.6]./5.2,linestyle=:dash,
      label = false)
=#
#=
f = h5open("EDcpgeL14OBCO0sigman2scaled.h5","r")
xs = read(f, "xs")
cpgeO0 = read(f, "cpge")
close(f)
f = h5open("EDcpgeL14OBCOn2sigman2.h5","r")
xs = read(f, "xs")
cpgeOn2 = read(f, "cpge")
close(f)
f = h5open("EDcpgeL14OBCOn3sigman2.h5","r")
xs = read(f, "xs")
cpgeOn3 = read(f, "cpge")
close(f)
f = h5open("EDcpgeL14OBCOn4sigman2.h5","r")
xs = read(f, "xs")
cpgeOn4 = read(f, "cpge")
close(f)
f = h5open("EDcpgeL14OBCO2n2sigman2.h5","r")
xs = read(f, "xs")
cpgeO2n2 = read(f, "cpge")
close(f)
f = h5open("EDcpgeL14OBCO3n2sigman2.h5","r")
xs = read(f, "xs")
cpgeO3n2 = read(f, "cpge")
close(f)
f = h5open("EDcpgeL14OBCO4n2sigman2.h5","r")
xs = read(f, "xs")
cpgeO4n2 = read(f, "cpge")
close(f)
f = h5open("EDcpgeL14OBCO5n2sigman2.h5","r")
xs = read(f, "xs")
cpgeO5n2 = read(f, "cpge")
close(f)
f = h5open("EDcpgeL14OBCO7n2sigman2.h5","r")
xs = read(f, "xs")
cpgeO7n2 = read(f, "cpge")
close(f)
f = h5open("EDcpgeL14OBCO9n2sigman2.h5","r")
xs = read(f, "xs")
cpgeO9n2 = read(f, "cpge")
close(f)

factor = 12*pi*0.01/5.28^4
plot(xs, cpgeOn2*factor, 
	linewidth = 2,
	xlabel = L"E_F/D",
     	ylabel = "CPGE"*L"*12\pi\Omega/D^4",
	#ylims = [-1,1],
	label=L"\Omega=10^{-2}",
	margin=3mm)
#plot!(xs, cpgeOn1, label=L"\Omega=10^{-1}")
plot!(xs, [cpgeO2n2*2*factor cpgeO3n2*3*factor cpgeO4n2*4*factor],
     label = [L"\Omega=2*10^{-2}" L"\Omega=3*10^{-2}" L"\Omega=4*10^{-2}"])
vline!([-1.6,-1.14,-1.124, -0.8,-0.352,-0.312,0,0.8,1.6]./5.28,linestyle=:dash,
      label = false,linecolor=:red)
vline!([-1.14-0.0528,-1.14-0.0528*3,-0.352-0.0528,-0.312+0.0528]./5.28,linestyle=:dash,
      label = false, linecolor=:black)
=#
#=
f = h5open("EDcpgeL14PBCOn2etan2.h5","r")
xs = read(f, "xs")
cpgeOn2 = read(f, "cpge")
y_all = read(f,"y_all")
close(f)
f = h5open("KPMcpgeL14NCs.h5","r")
NCs = read(f, "NCs")
Efs = read(f, "Efs")
cpgeKPMNCs = read(f, "CPGEs")
close(f)

sam = 4
maxED = maximum(cpgeOn2)
maxKPM = maximum(imag.(cpgeKPMNCs[:,sam]))
plot(xs, cpgeOn2/maxED,
        linewidth = 2,
        xlabel = L"E_F/D",
        ylabel = "CPGE scaled",
        label="ED",
	title = L"L=14",
        margin=3mm)
plot!(Efs, imag.(cpgeKPMNCs[:,:])/maxKPM,#label = L"sam="*"$sam")
     label = [L"N_C=16" L"N_C=32" L"N_C=64" L"N_C=128"])
vline!([-1.6,-1.14,-1.124, -0.8,-0.352,-0.312,0,0.8,1.6]./5.28,linestyle=:dash,
      label = false,linecolor=:red)
=#
#=
f = h5open("CPGEdata/E2s5L100M2W05gamma0Sam1kNoDis.h5","r")
E2s = read(f,"E2s")
#disorders = read(f,"disorders")
xs = read(f,"xs")
close(f)
scatter(xs, E2s[1:3,:]',label=false,
    xlabel = "Sample",
   ylabel = L"E^2",
   title = L"L=100,M=2,W=0.5",
   yaxis = :log,
   #ylims = [0,0.05], 
   framestyle = :box,grid=false,
    markerstrokewidth=0, markersize=3,
    xtickfontsize=12, ytickfontsize=12,
    xguidefontsize=12, yguidefontsize=12,
    legendfontsize=10,
    margin=2mm,
  )
#savefig("~/Desktop/E2sL100.pdf")
=#
#=
f = h5open("CPGEdata/ErareL18M2W05gamma04Sam1050TBCsPi.h5","r")
#f = h5open("CPGEdata/ErareL18M2W0gamma08TBCs.h5","r")
thetas = read(f, "thetas")
Esx = read(f, "Esx")
Esy = read(f, "Esy")
Esz = read(f, "Esz")
close(f)

l = @layout [a b]
p1 = plot(thetas, Esx[5770:5833,:]',label=false,
    xlabel = L"\theta_x",
   ylabel = L"E",
   #ylim = [-1,-0.6],
   title = "level 5770-5833", #5825-5840
  framestyle = :box,grid=false,
 xtickfontsize=12, ytickfontsize=12,
    xguidefontsize=12, yguidefontsize=12,
    legendfontsize=10,
    margin=2mm,
)
#p2 = plot(thetas, Esy[5825:5840,:]', label = false,
#    xlabel = L"\theta_y",
#   ylabel = L"E")
p3 = plot(thetas, Esz[5770:5833,:]',label = false,
    xlabel = L"\theta_z",
   ylabel = L"E",
     title = "level 5770-5833",
  framestyle = :box,grid=false,
 xtickfontsize=12, ytickfontsize=12,
    xguidefontsize=12, yguidefontsize=12,
    legendfontsize=10,
    margin=2mm,)
plot(p1,p3, layout = l)
=#
#=
f = h5open("EigenStatesRareM2g04TPBC.h5","r")
EsPBC = read(f, "EsPBC")
EsTBC = read(f, "EsTBC")
iprPBC = read(f, "iprPBC")
iprTBC = read(f, "iprTBC")
close(f)

xs = collect(1:11664)
scatter(xs, iprPBC,label="PBC",
    xlabel = "Eigen State",
    ylabel = "IPR"*L"(k)=\sum_i |\phi_k(i)|^4",
   title = L"L=18,M=2,\gamma=0.4,W=0.5",
   #ylims = [0,0.05],
   yaxis = :log,
   framestyle = :box,grid=false,
    markerstrokewidth=0.1, markersize=2.5,
    xtickfontsize=12, ytickfontsize=12,
    xguidefontsize=12, yguidefontsize=12,
    legendfontsize=10,
    legend = :topleft,
    margin=2mm,
  )
scatter!(xs, iprTBC,label="TBC",markerstrokewidth=0.1, markersize=2.5)
=#
