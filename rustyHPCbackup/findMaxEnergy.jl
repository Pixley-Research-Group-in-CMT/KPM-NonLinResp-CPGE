using Distributed
using ProgressMeter
using Random, KrylovKit #,Arpack
using SparseArrays
using JLD2,HDF5

f = h5open("/mnt/home/awu10/ceph/CPGEdata/E2s5L100M2W05gamma0Sam1k.h5","r")
E2s = read(f,"E2s")
disorders = read(f,"disorders")
xs = read(f,"xs")
samples = read(f,"samples")
close(f)

fn = "/mnt/home/awu10/ceph/CPGEdata/spH0L100M2gam0.jld2"
H0 = load(fn,"H0")
Ls = load(fn,"Ls")
L = Ls[1]
M = load(fn,"M")
γ = load(fn,"gamma")
W = 0.5

Emax = zeros(samples)
Threads.@threads for k in 1:samples
        Ham = H0 + spdiagm(disorders[:,k])
        #@time Etemp, Vtemp = eigs(Ham;nev=1, which=:SM,tol=1e-3)
        @time Etemp, Vtemp = eigsolve(Ham,2) #(Ham;nev=1, which=:SM,tol=1e-3,maxiter=10)
        Emax[k] = Etemp[1] #real.(eigvals(Ham*Ham))
        println(k)
        #sleep(0.01)
end

f = h5open("EmaxL100M2W05gamma0Sam1k.h5","w")
write(f, "Emax",Emax)
close(f)



