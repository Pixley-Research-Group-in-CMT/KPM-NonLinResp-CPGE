using Distributed
#using ProgressMeter
using Random, KrylovKit #,Arpack
using JLD2,HDF5
using Measures

include("model.jl")

@show BLAS.get_num_threads()
#BLAS.set_num_threads(1)
#@show BLAS.get_num_threads()
#=
L = 100
Ls = [L,L,L]
M = 2.0 # 1.5 2.0
t = 1.0
γ = 0.1  # 0.8
W = 0.5

twists = [pi,pi,pi]
boundary = exp.(1im.*twists)
samples = 350
H0 = TwoBandModel(Ls, t, M, γ;boundary=boundary)

disorders = zeros(2*L^3,samples)
E2s = zeros(2*L^3, samples)
Threads.@threads for k in 1:samples
	disorders[:, k] = generateDisorder(Ls, W)
	Ham = Matrix(H0 + spdiagm(disorders[:,k]))
	E2s[:,k] = eigvals(Ham*Ham)
	println(k)
	#sleep(0.01)
end

#using Plots
xs = collect(1:350)
#scatter(xs, E2s[1:10,:]')

f = h5open("/mnt/home/awu10/ceph/CPGEdata/E2sL18M15W05gamma0Sam1050.h5","r")
E2s = read(f,"E2s")
disorders = read(f,"disorders")
close(f)


# find rare state
disind = argmin(E2s)[2]
thetas = zeros(100)
Esx = zeros(2*L^3, 100)
Esy = zeros(2*L^3, 100)
Esz = zeros(2*L^3, 100)
# Turn off disorder
#dimx, dimy = size(disorders)[1], size(disorders)[2]
#disorders = zeros(dimx,dimy)
Threads.@threads for k in 1:100
	thetas[k] = k*2*pi/100
	twist = [thetas[k], 0, 0]
	boundary = exp.(1im.*twist)
	H1 = TwoBandModel(Ls, t, M, γ;boundary=boundary)
	Ham = Matrix(H1 + spdiagm(disorders[:,disind]))
	Esx[:,k] = eigvals(Ham)
	
	twist = [0,thetas[k],0]
        boundary = exp.(1im.*twist)
	H1 = TwoBandModel(Ls, t, M, γ;boundary=boundary)
        Ham = Matrix(H1 + spdiagm(disorders[:,disind]))
        Esy[:,k] = eigvals(Ham)
	
	twist = [0, 0, thetas[k]]
        boundary = exp.(1im.*twist)
        H1 = TwoBandModel(Ls, t, M, γ;boundary=boundary)
        Ham = Matrix(H1 + spdiagm(disorders[:,disind]))
        Esz[:,k] = eigvals(Ham)

	println(k)
	#sleep(0.01)
end

#f = h5open("ErareL18M2W05gamma08Sam1050TBCsPi.h5","w")
f = h5open("ErareL18M2W05gamma08Sam1050TBCs.h5","w")
write(f,"Esx",Esx)
write(f,"Esy",Esy)
write(f,"Esz",Esz)

close(f)
=#
#=
pt = "/home/aw666/data/CPGEdata/"
fn = pt*"E2s5L100M2W05gamma0Sam100TBC.jld2"
E2s = load(fn,"E2s")
disind = argmin(E2s[1,:]) #argmax(E2s[1,:])#argmin(E2s[1,:])
disorder = load(fn,"disorders")[:,disind]
Nvs = 5
E2sxyz = zeros(Nvs,3)
IPRsxyz = zeros(Nvs, 3)

tind = parse(Int64,ARGS[1])
theta = tind*2*pi/100

function getEsIPRs(twist,disorder)
	boundary = exp.(1im.*twist)
	H0 = TwoBandModel(Ls, t, M, γ;boundary=boundary)
	Ham = H0 + spdiagm(disorder)
	Ham = Ham * Ham
	@time Etemp, Vtemp,_ = eigsolve(Ham,Nvs,:SR)
	E2s = abs.(Etemp[1:Nvs])
	Ipr = zeros(Nvs)
	for ind = 1:Nvs
		Ipr[ind] = sum((abs.(Vtemp[ind])).^4)
	end
	return E2s, Ipr
end

E2sxyz[:,1], IPRsxyz[:,1] = getEsIPRs([theta,0,0],disorder)
E2sxyz[:,2], IPRsxyz[:,2] = getEsIPRs([0,theta,0],disorder)
E2sxyz[:,3], IPRsxyz[:,3] = getEsIPRs([0,0,theta],disorder)

pt = "/scratch/aw666/CPGEdata/"  # amarel
fn = pt*"E2sIPRrareL100M2W05gamma01Twist"*ARGS[1]*".jld2"
#fn = pt*"E2sIPRpertL100M2W05gamma04Twist"*ARGS[1]*".jld2"
jldsave(fn; disind=disind,E2sxyz=E2sxyz,IPRsxyz=IPRsxyz)
=#

t = 1
M = 2
L = 18
Ls = [L,L,L]
fn = "CPGEdata/E2s5L18M2W05gamma0Sam1000.jld2"
#f = h5open(fn,"r")
E2s = load(fn,"E2s") #read(f,"E2s")
RareInd = argmin(E2s[1,:])
PertInd = argmax(E2s[1,:])
disorders = load(fn,"disorders") #read(f,"disorders")
raredis = disorders[:,RareInd]
pertdis = disorders[:,PertInd]
γs = [0,0.1,0.2,0.3,0.4,0.5]
Es = zeros(2*L^3,6)
Vs = zeros(ComplexF64,2*L^3,2*L^3,6)
IPRs = zeros(2*L^3,6)
println("start to compute")
for k = 1:6
	boundary = [-1,-1,-1]
	H0 = TwoBandModel(Ls, t, M, γs[k];boundary=boundary)
        disorder = raredis
	Ham = H0 + spdiagm(disorder)
	@time Etemp, Vtemp = eigen(Matrix(Ham))
	Es[:,k] = Etemp
	Vs[:,:,k] = Vtemp
	for ind = 1:(2*L^3)
                IPRs[ind,k] = sum((abs.(Vtemp[:,ind])).^4)
        end
	println(k)
end

function indToPos(ind,L)
	ind = div(ind+1,2)
	rz = (ind-1) % L + 1
c	ind = div((ind - rz),L)
	ry = (ind-1)%L + 2
	ind = div((ind - ry),L)
	rx = (ind-1)%L + 2
	return [rx,ry,rz]
end

function WFVSr(psi)
	center = indToPos(argmax(abs.(psi)),L)
	distmaps = Dict()
	for k in 1:2*L^3
		vector = indToPos(k,L).-center
		dist = sum(vector.^2)
		if haskey(distmaps,dist)
			push!(distmaps[dist],k)
		else
			distmaps[dist] = [k]
		end
	end
	rs = []
	ys = []
	for (key,value) in distmaps
		push!(rs,sqrt(key))
		push!(ys,maximum(abs.(psi[value])))
	end
	inds = sortperm(rs)
	rs = rs[inds]
	ys = ys[inds]
	rs = convert(Array{Float64,1}, rs)
	ys = convert(Array{Float64,1}, ys)
	return rs, ys
end
ind = 5833
psi = Vs[:,ind,1]
rs, ys = WFVSr(psi)
endind = 100 #length(rs)
xs = log.(rs[2:endind])
Xs = ones(endind-1,2)
Xs[:,1] = xs
y = log.(ys[2:endind])
b = Xs \ y
ypred = Xs*b
plot(exp.(xs),exp.([y ypred]),xaxis=:log,yaxis=:log)

rss = []
yss = []
for k = 1:6
	psitemp = Vs[:,ind,k]
	rt, yt = WFVSr(psitemp)
	push!(rss,rt[2:length(rt)])
	push!(yss,yt[2:length(rt)])
end

plot(exp.(xs),exp.(ypred),linewidth=4,color=:red,label=false,
     xaxis=:log,yaxis=:log,
     ylabel=L"\sqrt{\psi^\dagger\psi}",
     xlabel=L"r",
	xlim=(0.8,12),
	xticks=([10^0,10^0.5,10^1],
		[L"10^0",L"10^{0.5}",L"10^1",L"1333"]),
	yticks=([10^(-3),10^(-2),10^(-1)],
		[L"10^{-3}",L"10^{-2}",L"10^{-1}"]),
	ylim = (10^(-3),0.3),
	legend =:bottomleft,
        framestyle = :box,grid=false,
        xtickfontsize=12, ytickfontsize=12,
        xguidefontsize=12, yguidefontsize=12,
        legendfontsize=12,margin=5mm
    )
scatter!(rss[1],yss[1],markerstrokewidth=0, markersize=4,label=L"\gamma=0")
scatter!(rss[2],yss[2],markerstrokewidth=0, markersize=2,label=L"\gamma=0.1")
scatter!(rss[3],yss[3],markerstrokewidth=0, markersize=2,label=L"\gamma=0.2")
scatter!(rss[4],yss[4],markerstrokewidth=0, markersize=2,label=L"\gamma=0.3")
scatter!(rss[5],yss[5],markerstrokewidth=0, markersize=2,label=L"\gamma=0.4")
scatter!(rss[6],yss[6],markerstrokewidth=0, markersize=2,label=L"\gamma=0.5")
#display(Plots.plot!())

#savefig("~/Desktop/PsiVsRgammas.pdf")


scatter(E2s[1,:],
        label = false,
        xlabel=L"\mathrm{disorder\;realizations}",
        ylabel=L"E_n^2\;(\!\times 10^{-3})",
        legend =false,
        framestyle = :box,grid=false,
        xtickfontsize=12, ytickfontsize=12,
        xguidefontsize=14, yguidefontsize=14,
        legendfontsize=12,margin=5mm
        )
vline!([argmin(E2s[1,:]) argmax(E2s[1,:])],linewidth=1,linestyle=:dash,color=[:red :black],label=false)

scatter(IPRs[:,[1,4,5]],
        label = [L"\gamma=0" L"\gamma=0.3" L"\gamma=0.4"],
	xlabel=L"\mathrm{state}\; n",
        ylabel=L"\mathrm{IPR}(n)=\sum_i |\psi_n(i)|^4",
	markerstrokewidth=0, markersize=[4 3 3],
        yaxis=:log,
        legend =:topright,
        framestyle = :box,grid=false,
        xtickfontsize=12, ytickfontsize=12,
        xguidefontsize=14, yguidefontsize=14,
        legendfontsize=12,
        )
vline!([5833],linewidth=1,linestyle=:dash,color=:black,label=false)



