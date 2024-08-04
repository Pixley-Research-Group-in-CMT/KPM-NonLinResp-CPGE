using Statistics
using LinearAlgebra, Arpack, SparseArrays
using KPM
using JLD2#,HDF5
using CUDA

include("model.jl")
@show CUDA.has_cuda()
@show KPM.whichcore()
# check number of threads
println("Number of threads: $(Threads.nthreads())")

L = 18 #100
Ls = [L,L,L]
M = 2.0 # 2.0
t = 1.0
γ = 0.3  # 0.8
W = 0.0 #sqrt(0.01)
twists = 2*pi*rand(3) #[pi,pi,pi]
@show twists
boundary = exp.(1im.*twists)
@time H0 = TwoBandModel(Ls, t, M, γ;boundary=boundary)

disind = parse(Int64,ARGS[1])
#disorder = generateDisorder(Ls, W)
#disorder = generatePotDisorder(Ls, W)

pt = "/home/aw666/data/CPGEdata/"
fn = pt*"E2s5L18M2W05gamma0Sam1000.jld2"
E2s = load(fn,"E2s")
disk = disind #argmax(E2s[1,:])#argmin(E2s[1,:])
disorder = load(fn,"disorders")[:,disk]

#disorder = zeros(2*L^3)
Ham = H0 + spdiagm(disorder) #spdiagm(disorders[:,disind])
Jx, Jy, Jz = velocityOperator(Ls,Ham) # same as no disorder

D = 6.0 #(M+2)+sqrt(1+gamma^2)
H_norm = Ham./D

pt = "/scratch/aw666/CPGEdata/"  # amarel
fn = pt*"HdisL$(L)gamma03W05PotDisTBC"*"$disind"*".jld2"
jldsave(fn; H_norm=H_norm,D=D, Jx=Jx,Jy=Jy, Jz=Jz)

