using Distributed
#using ProgressMeter
using Random, KrylovKit #,Arpack
using JLD2 #,HDF5

include("model.jl")

@show Threads.nthreads()
#@show BLAS.get_num_threads()
#BLAS.set_num_threads(1)
#@show BLAS.get_num_threads()

L = 18
Ls = [L,L,L]
M = 2.0 # 2.0
t = 1.0
gamma = 0.0  # 0.8
twists = [pi,pi,pi]
boundary = exp.(1im.*twists)
@time H0 = TwoBandModel(Ls, t, M, gamma;boundary=boundary)

#=
fn = "/mnt/home/awu10/ceph/CPGEdata/spH0L100M2gam0.jld2"
H0 = load(fn,"H0")
Ls = load(fn,"Ls")
L = Ls[1]
M = load(fn,"M")
γ = load(fn,"gamma")
=#
samples = 1000
W = 0.5
Nvs = 5

disorders = zeros(2*L^3,samples)
E2s = zeros(Nvs, samples)
IPRs = zeros(Nvs, samples)
Threads.@threads for k in 1:samples
#for k in 1:samples
	disorders[:, k] = generatePotDisorder(Ls, W)
        #Ham = Matrix(H0 + spdiagm(disorders[:,k]))
	Ham = H0 + spdiagm(disorders[:,k])
	#@time Etemp, Vtemp = eigs(Ham;nev=1, which=:SM,tol=1e-3)
	Ham = Ham * Ham
	@time Etemp, Vtemp,_ = eigsolve(Ham,5,:SR) #(Ham;nev=1, which=:SM,tol=1e-3,maxiter=10)
	E2s[:,k] = abs.(Etemp[1:Nvs]) #real.(eigvals(Ham*Ham))
	for ind = 1:Nvs
		IPRs[ind,k] = sum((abs.(Vtemp[ind])).^4)
	end
	println(k)
        #sleep(0.01)
end

xs = collect(1:samples)
#pt = "/scratch/aw666/CPGEdata/"  # amarel
pt = "/home/aw666/data/CPGEdata/"
#fn = pt*"E2s5L$(L)M2W05gamma0Sam$(samples)Ind"*ARGS[1]*"MyTBC.jld2"
fn = pt*"E2s5L$(L)M2W05gamma0Sam$(samples).jld2"
jldsave(fn; disorders=disorders,E2s=E2s,IPRs=IPRs,xs=xs)

#=
fn = "E2s5L100M2W05gamma0Sam$(samples)"*ARGS[1]*"MyTBC.h5"
f = h5open(fn,"w")
write(f, "disorders", disorders)
write(f, "E2s",E2s)
write(f,"samples", samples)
write(f,"xs", xs)
close(f)
=#




