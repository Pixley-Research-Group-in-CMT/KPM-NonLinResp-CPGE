using Plots
using LaTeXStrings
using Measures
using KPM
using JLD2

include("model.jl")

M = 2.0 # 2.0
t = 1.0
γ = 0.8 # 0.8
Npts = 1000
fn = "CPGEdata/DOSmuL100gam08D6NC215NR10TBC100.jld2"
mu = load(fn,"mu")
D = load(fn,"D")
NC = 2^15

kz = collect(LinRange(-0.5,0.5,Npts))
Ekp,Ekn = zeros(Npts), zeros(Npts)
errp,errn = zeros(Npts), zeros(Npts)
for k = 1:Npts
	Ep, En = dispersion(0, 0,kz[k], t,M,γ)
	Ekp[k] = Ep
	Ekn[k] = En
	errp[k] = KPM._dos_single(mu, D, Ep, NC)
	errn[k] = KPM._dos_single(mu, D, En, NC)
end
cors = palette(:default)
plot(pi*kz,Ekp,ribbon=(errp,errp),color=cors[1],linewidth=2, fillalpha=0.3)
plot!(pi*kz,Ekn,ribbon=(errn, errn),color=cors[1],linewidth=2, fillalpha=0.3)
hline!([0.6], color=:black, 
       linestyle=:dot, linewidth=2, 
       fill=(-1.45, "#8ACB88"), fillalpha=0.3,
       xticks = [], yticks=[],grid=false,
      showaxis=false, legend=false)
#savefig("~/Desktop/Fig1Dispersion.pdf")

