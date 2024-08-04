using Distributed
using ProgressMeter
using Random
using HDF5

include("model.jl")

@show BLAS.get_num_threads()
#BLAS.set_num_threads(1)
#@show BLAS.get_num_threads()

L = 18
Ls = [L,L,L]
M = 2.0 # 1.5 2.0
t = 1.0
γ = 0.8  # 0.8
W = 0.5

#=
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
=#

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
write(f,"disind", disind)
write(f,"thetas", thetas)
close(f)



