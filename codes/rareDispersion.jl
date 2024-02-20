using Plots
using LaTeXStrings
using Measures

include("model.jl")

function RhoAdispersion(γ)
	t = 1
	M = 2
	Nsam = 1000
	kz1s = collect(LinRange(-0.5,0.5,Nsam))
	Eskz = zeros(Nsam,2)
	for k = 1:Nsam
		E1, E2 = dispersion(0,0, kz1s[k], t, M, γ)
		Eskz[k,:] = [E1,E2]
	end
	Nsam = 100
	kxs = collect(LinRange(-0.5,0.5,Nsam))
	kys = copy(kxs)
	kzs = copy(kxs)
	Etot = zeros(Nsam,Nsam,Nsam,2)
	for x=1:Nsam, y=1:Nsam, z=1:Nsam
		E1, E2 = dispersion(kxs[x],kys[y], kzs[z], t, M, γ)
		Etot[x,y,z,:] = [E1,E2] 
	end
	Etot = reshape(Etot,1,:)
	Etot = vec(Etot)
	sort!(Etot)
	return kz1s, Eskz, Etot
end
γ = 0.0
kz1s, Eskz0,Etot0 = RhoAdispersion(γ)
γ = 0.4
kz1s, Eskz02,Etot02 = RhoAdispersion(γ)


cors = palette(:default)
histogram([Etot0 Etot02],bins=1000,normalize=:pdf,
	  xlabel = L"E",
	  ylabel = L"\rho(E)",
	  xlim=(-1.5,1.5),
	  yaxis=:log,ylim=(10^(-4.2),1),
	  label = [L"\gamma=0" L"\gamma=0.4"],
	  legend = :bottomleft,
     	framestyle = :box,grid=false,
 	xtickfontsize=12, ytickfontsize=12,
    	xguidefontsize=12, yguidefontsize=12,
    	legendfontsize=10,
    	margin=3mm
	  )
vline!([-0.4 0 0.4],linewidth=1,linestyle=:dash,color=:black,label=false)
plot!(kz1s.*2, [Eskz0 Eskz02],
      xlabel = L"k_z/\pi",
      ylabel = L"E(k_z)",
      color = [cors[1] cors[1] cors[2] cors[2]],
        inset = (1, bbox(0.68,0.5,0.3,0.37)), subplot = 2,
        legend=false,
        framestyle = :box,grid=false,
        xtickfontsize=10, ytickfontsize=10,
        xguidefontsize=10, yguidefontsize=10,
    )
hline!([-0.4 0 0.4],linewidth=1,linestyle=:dash,color=:black,label=false,
      subplot = 2,)
#savefig("~/Desktop/rhoDispersionG0G02.pdf")
