using Distributed
using ProgressMeter
using Random, KrylovKit #,Arpack
using JLD2,HDF5

include("model.jl")

#@show BLAS.get_num_threads()
#BLAS.set_num_threads(1)
#@show BLAS.get_num_threads()

L = 100
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
samples = 200
W = 0.5

disorders = zeros(2*L^3,samples)
E2s = zeros(5, samples)
Threads.@threads for k in 1:samples
        disorders[:, k] = generateDisorder(Ls, W)
        #Ham = Matrix(H0 + spdiagm(disorders[:,k]))
	Ham = H0 + spdiagm(disorders[:,k])
	#@time Etemp, Vtemp = eigs(Ham;nev=1, which=:SM,tol=1e-3)
	Ham = Ham * Ham
	@time Etemp, Vtemp = eigsolve(Ham,5,:SR) #(Ham;nev=1, which=:SM,tol=1e-3,maxiter=10)
	E2s[:,k] = Etemp[1:5] #real.(eigvals(Ham*Ham))
        println(k)
        #sleep(0.01)
end

xs = collect(1:samples)

fn = "E2s5L100M2W05gamma0Sam200"*ARGS[1]*"MyTBC.h5"

f = h5open(fn,"w")
write(f, "disorders", disorders)
write(f, "E2s",E2s)
write(f,"samples", samples)
write(f,"xs", xs)
close(f)




