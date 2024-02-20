using Plots
using LaTeXStrings
using Measures
using JLD2
using Statistics
#=
histogram(Eana,bins=1000,normalize=:pdf,
	  label = "analytical "*L"E",
     xlabel = L"E",
     ylabel = L"\rho(E)",
     title = L"L=200,\gamma=0.4"*", PBC",
     legend = :bottomright,
     framestyle = :box,grid=false,
 xtickfontsize=12, ytickfontsize=12,
    xguidefontsize=12, yguidefontsize=12,
    legendfontsize=10,
    margin=3mm,)
histogram!(Eana,bins=1000,normalize=:pdf,
	   xlim = [-1.2,1.2],
    framestyle = :box,
    inset = (1, bbox(0.1,0.5,0.35,0.4)), subplot = 2,
    legend = false,
    xtickfontsize=10, ytickfontsize=10,
    xguidefontsize=10, yguidefontsize=10,grid=false
)
vline!([-0.4,0.4], linestyle=:dash,subplot=2)
savefig("~/Desktop/DOSL200gam04.pdf")
=#
#=
pt = "./CPGEdata/"
fn = pt*"GammaLsNC64CPUrunsams100.jld2"
runtimesCPU = load(fn, "runtimes")
fn = pt*"GammaLsNC64GPUrunsams100.jld2"
runtimesGPU = load(fn, "runtimes")
Ls = load(fn,"Ls")
CPUmedian = vec(median(runtimesCPU, dims=1))
GPUmedian = vec(median(runtimesGPU, dims=1))
CPUmean = vec(mean(runtimesCPU, dims=1))
GPUmean = vec(mean(runtimesGPU, dims=1))
cors = palette(:default)
plot(Ls, [CPUmedian CPUmean GPUmedian GPUmean],
     markershape=[:circle :utriangle],
     color=[cors[1] cors[1] cors[2] cors[2]],
     linestyle = [:solid :dash :solid :dash],
    yaxis=:log,
    xlabel = L"L",
    ylabel = L"\mathrm{run\; time}(s)",
    labels = [L"\mathrm{CPU\; median}" L"\mathrm{CPU\;mean}" L"\mathrm{GPU\;median}" L"\mathrm{GPU\;mean}"],
    title=L"N_R=1,N_C=2^6",
        legend = :topleft,
        framestyle = :box,grid=false,
        xtickfontsize=12, ytickfontsize=12,
        xguidefontsize=12, yguidefontsize=12,
        legendfontsize=12,
        margin=2mm,
   )
cors = palette(:grays,10)
plot!(Ls,[CPUmedian./GPUmedian CPUmean./GPUmean],
      color=[cors[1] cors[8]],markershape=:circle,
    labels = [L"\mathrm{median}" L"\mathrm{mean}"],  
    ylabel=L"\mathrm{speedup}",
        inset = (1, bbox(0.5,0.3,0.45,0.4)), subplot = 2,
     legend =:topleft,
     framestyle = :box,grid=false,
 xtickfontsize=10, ytickfontsize=10,
    xguidefontsize=10, yguidefontsize=10,
    legendfontsize=12,
    )
#savefig("~/Desktop/GPUCPUcpge.pdf")
=#
#=
pt = "./CPGEdata/"
fn = pt*"E2s5L100M2W05gamma0Sam100TBCnoDis.jld2"
E2s = load(fn, "E2s")
IPRs = load(fn, "IPRs")
scatter(1000*E2s[1:2,:]',
	label = [L"E_1" L"E_2"],
	xlabel=L"\mathrm{disorder\;realizations}",
	ylabel=L"E_n^2\;(\!\times 10^{-3})",
	#title=L"L=100,W=0.5",
	ylim = [0,1.75],
	yticks = ([0.0,0.5,1.0,1.5],
                       [L"0.0",L"0.5", L"1.0",L"1.5"]),
	legend =:bottomleft,
     	framestyle = :box,grid=false,
	xtickfontsize=12, ytickfontsize=12,
	xguidefontsize=14, yguidefontsize=14,
	legendfontsize=12,
	)
vline!([argmin(E2s[1,:]) argmax(E2s[1,:])],linewidth=1,linestyle=:dash,color=[:red :black],label=false)
#savefig("~/Desktop/E2sAlldis.pdf")
scatter(IPRs[1:2,:]',
        label = [L"E_1" L"E_2"],
        xlabel=L"\mathrm{disorder\;realizations}",
	ylabel=L"\mathrm{IPR}(n)=\sum_i |\psi_n(i)|^4",
        yaxis=:log,
	legend =:right,
        framestyle = :box,grid=false,
        xtickfontsize=12, ytickfontsize=12,
        xguidefontsize=14, yguidefontsize=14,
        legendfontsize=12,
        )
vline!([argmin(E2s[1,:]) argmax(E2s[1,:])],linewidth=1,linestyle=:dash,color=[:red :black],label=false)
#savefig("~/Desktop/IPRsAlldis.pdf")
=#
#=
pt = "./CPGEdata/"
fn = pt*"E2sIPRpertL100M2W05gamma0Twists100.jld2"
E2sxpert = load(fn,"E2sx")
fn = pt*"E2sIPRrareL100M2W05gamma0Twists100.jld2"
E2sxrare = load(fn,"E2sx")
thetas = collect(1:100)*2*pi/100
plot(thetas, 1000*E2sxrare[1:2,:]',
     label = [L"E_1" L"E_2"],
     xlabel=L"\theta_x",
        ylabel=L"E^2_n",
        #yaxis=:log,
        legend =:topleft,
        framestyle = :box,grid=false,
        xtickfontsize=12, ytickfontsize=12,
        xguidefontsize=14, yguidefontsize=14,
        legendfontsize=12,size=(500, 400)
    )
plot!(thetas, 1000*E2sxpert[1:2,:]',
      title=L"\mathrm{perturbative}",
	inset = (1, bbox(0.65,0.1,0.35,0.3)), subplot = 2,
        legend=false,
	framestyle = :box,grid=false,    
        xtickfontsize=12, ytickfontsize=12,
        xguidefontsize=14, yguidefontsize=14,
        legendfontsize=12,
    )
=#
