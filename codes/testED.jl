using Distributed
using ProgressMeter
using Random
using HDF5
#using Plots
#pyplot()

include("model.jl")
include("EDcpge.jl")
#=
function IPR(Vs)
        n = size(Vs)[1]
        res = zeros(n,1)
        for k in 1:n
                res[k] = sum((abs.(Vs[:,k])).^4)
        end
        return res
end
=#

L = 200
Ls = [L,L,L]
M = 2.0 # 2.0
t = 1.0
γ = 0.8 # 0.8
@show Emax = (M+2)+sqrt(1+γ^2)
#Eana = spectrum(Ls,t,M,γ)

#Emax *= 1.5 # extend beyond Emax 

pts = 1000
paths = [[0.5,0.5,0.5],[0.0,0.0,0.0],[0.5,-0.5,0.0],[0.375,-0.375,0.0],[0.0,0.0,0.0]]
paths = [[0.0, 0.0, -0.5],[0.0,0.0,0.5]]
pathL = size(paths)[1]
kpoints = zeros(3,pts*(pathL-1))
for k in 1:pathL-1
	kpoints[:,(pts*(k-1)+1):(pts*k)] =  pathInMomentum(paths[k],paths[k+1];points = pts)
end
Eks = zeros(pts*(pathL-1),2)
for k in 1:(pts*(pathL-1))
	Ep,En = dispersion(kpoints[1,k], kpoints[2,k],kpoints[3,k], t,M,γ)
	Eks[k, :] = [Ep, En]
	if abs(Ep-En-0.8) < 0.005
		println(Ep," ", En, " ",kpoints[3,k])
	end
end

#=
boundary = [1,1,1] #exp.(1im*2*pi*rand(3))
H0 = TwoBandModel(Ls,t,M,γ;boundary=boundary)
H = H0./Emax
#@time Es, Vs = eigen(Matrix(H0))
#Eana = spectrum(Ls,t,M,γ)

hx, hy, hz = velocityOperator(Ls,H0)

@time Es, hxp, hyp, hzp = EDvelocity(Matrix(H), hx, hy, hz)
h = 0.002
x_all = collect(-1:h:1)
y_all = complex(x_all)

Ω = 2*10^(-2)
w1 = 0.8/Emax
@show w2 = -w1
@show w1 += Ω
@show σ = 10^(-2)
η = σ
@show ef = -1.143
@time v1 = EDcpge(ef/Emax, w1, w2, Es, hxp, hyp, hzp;σ=σ,η=η)
@time v2 = EDcpge(ef/Emax, w2, w1, Es, hyp, hxp, hzp;σ=σ,η=η)
=#
#=
@showprogress 1 "Computing" for (i, x) in enumerate(x_all)
	y_all[i] = EDcpge(x, w1, w2, Es, hxp, hyp, hzp;σ=σ,η=η)+EDcpge(x, w2, w1, Es, hyp, hxp, hzp;σ=σ, η=η)
	#println(i)
	sleep(0.001)
end
@show maximum(abs.(y_all))
y_int = cumsum(y_all)*h
maximum(abs.(y_int))
=#
#=
f = h5open("E2sL18M15W05gamma0Sam1050.h5","r")
E2s = read(f,"E2s")
disorders = read(f,"disorders")
close(f)

disind = argmin(E2s)[2]
boundary = [1,1,1]
H1 = TwoBandModel(Ls, t, M, γ;boundary=boundary)
Ham = Matrix(H1 + spdiagm(disorders[:,disind]))
@time EsPBC, VsPBC = eigen(Ham)
@time iprPBC = IPR(VsPBC)

boundary = [-1,-1,-1]
H1 = TwoBandModel(Ls, t, M, γ;boundary=boundary)
Ham = Matrix(H1 + spdiagm(disorders[:,disind]))
@time EsTBC, VsTBC = eigen(Ham)
@time iprTBC = IPR(VsTBC)

f = h5open("EigenStatesRareTPBC.h5","w")
write(f,"EsPBC",EsPBC)
write(f, "VsPBC",VsPBC)
write(f, "EsTBC",EsTBC)
write(f,"VsTBC",VsTBC)
write(f,"iprPBC",iprPBC)
write(f,"iprTBC", iprTBC)
write(f,"thetas", thetas)
close(f)
=#
#p1 = scatter(iprPBC, yaxis=:log,markerstrokewidth=0.1,markersize=3)
#p2 = scatter(iprTBC, yaxis=:log,markerstrokewidth=0.1,markersize=3,reuse=false)
#=
xs = collect(1:11664)
scatter(xs, iprPBC,label="PBC",
    xlabel = "Sample",
    ylabel = "IPR"*L"(k)=\sum_i |\phi_k(i)|^4",
   title = L"L=18,M=1.5,W=0.5",
   #ylims = [0,0.05],
   yaxis = :log,
   framestyle = :box,grid=false,
    markerstrokewidth=0.1, markersize=3,
    xtickfontsize=12, ytickfontsize=12,
    xguidefontsize=12, yguidefontsize=12,
    legendfontsize=10,
    legend = :topleft,
    margin=2mm,
  )
scatter!(xs, iprTBC,label="TBC",markerstrokewidth=0.1, markersize=3)
=#
