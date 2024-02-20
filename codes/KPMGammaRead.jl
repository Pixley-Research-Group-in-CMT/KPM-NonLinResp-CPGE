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

f = h5open("/mnt/home/awu10/ceph/CPGEdata/E2s5L100M2W05gamma0Sam1k.h5","r")
E2s = read(f,"E2s")
disorders = read(f,"disorders")
xs = read(f,"xs")
samples = read(f,"samples")
close(f)
#=
fn = "/mnt/home/awu10/ceph/CPGEdata/spH0L100M2gam08PBC.jld2"
H0 = load(fn,"H0")
Jx = load(fn,"Jx")
Jy = load(fn,"Jy")
Jz = load(fn,"Jz")
Ls = load(fn,"Ls")
L = Ls[1]
M = load(fn,"M")
γ = load(fn,"gamma")
=#

L = 100
Ls = [L,L,L]
M = 2.0 # 2.0
t = 1.0
γ = 0.8  # 0.8
W = 0.5
twists = 2*pi*rand(3) #[pi,pi,pi]
@show twists
boundary = exp.(1im.*twists)
@time H0 = TwoBandModel(Ls, t, M, γ;boundary=boundary)

NC = 512
NR = 1

#disind = argmin(E2s[1,:])
#disind = argmax(E2s[1,:])
disind = parse(Int64,ARGS[1])
disorder = generateDisorder(Ls, W)
Ham = H0 + spdiagm(disorder) #spdiagm(disorders[:,disind])
Jx, Jy, Jz = velocityOperator(Ls,Ham) # same as no disorder

D = 6.0 #(M+2)+sqrt(1+gamma^2)
H_norm = Ham./D


NH = H_norm.n
mu_3d_xyz = zeros(ComplexF64, NC, NC, NC)
psi_in_l = exp.(2pi * 1im * rand(H_norm.n, NR));
KPM.normalize_by_col(psi_in_l, NR)
psi_in_r = psi_in_l
KPM.kpm_3d!(H_norm, Jx, Jy, Jz, NC, NR, NH, mu_3d_xyz, psi_in_l, psi_in_r; verbose=1, arr_size=16);


mu_3d_yxz = zeros(ComplexF64, NC, NC, NC)
psi_in_l = exp.(2pi * 1im * rand(H_norm.n, NR));
KPM.normalize_by_col(psi_in_l, NR)
psi_in_r = psi_in_l
KPM.kpm_3d!(H_norm, Jy, Jx, Jz, NC, NR, NH, mu_3d_yxz, psi_in_l, psi_in_r; verbose=1, arr_size=16);

#fn = "GammamnpL100gam02D6NC512W05maxNRs5.jld2"
fn = "/mnt/home/awu10/ceph/CPGEdata/GammamnpL100gam08D6NC512W05maxNR1Sam"*"$disind"*".jld2"
jldsave(fn; Ls=Ls, NC=NC, D = D, Gammaxyz=mu_3d_xyz, Gammayxz=mu_3d_yxz)





