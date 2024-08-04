using Statistics
using LinearAlgebra, Arpack, SparseArrays
using KPM
using JLD2,HDF5
using CUDA

include("model.jl")
@show CUDA.has_cuda()
@show KPM.whichcore()
# check number of threads
println("Number of threads: $(Threads.nthreads())")

L = 100
Ls = [L,L,L]
M = 2.0 # 2.0
t = 1.0
γ = 0.8  # 0.8
W = sqrt(0.35)
twists = 2*pi*rand(3) #[pi,pi,pi]
@show twists
boundary = exp.(1im.*twists)
@time H0 = TwoBandModel(Ls, t, M, γ;boundary=boundary)

disind = parse(Int64,ARGS[1])
disorder = generateDisorder(Ls, W)
Ham = H0 + spdiagm(disorder) #spdiagm(disorders[:,disind])
Jx, Jy, Jz = velocityOperator(Ls,Ham) # same as no disorder

D = 6.0 #(M+2)+sqrt(1+gamma^2)
H_norm = Ham./D

fn = "HdisL100gam08W2035dis"*"$disind"*".jld2"
jldsave(fn; H_norm=H_norm, Jx=Jx,Jy=Jy, Jz=Jz)

