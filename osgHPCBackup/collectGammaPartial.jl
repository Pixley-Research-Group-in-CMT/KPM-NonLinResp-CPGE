using Statistics
using LinearAlgebra, Arpack, SparseArrays
using KPM
using JLD2
using CUDA

@show CUDA.has_cuda()
@show KPM.whichcore()
# check number of threads
println("Number of threads: $(Threads.nthreads())")

disind = parse(Int64,ARGS[1])
#pt = "/public/angkun/CPGEdata/"  # osg
#fn = pt*"HdisL100gam08W2002dis"*"$disind"*".jld2"
fn = "HdisL100gam08W2035dis"*"$disind"*".jld2"
H_norm = load(fn, "H_norm")
Jx = load(fn, "Jx")
Jy = load(fn, "Jy")
Jz = load(fn, "Jz")

NC = 512
NR = 1
NH = H_norm.n
D = 6.0 #(M+2)+sqrt(1+gamma^2)
L = 100
Ls = [L,L,L]

mu_3d_xyz = zeros(ComplexF64, NC, NC, NC)
psi_in_l = exp.(2pi * 1im * rand(H_norm.n, NR));
KPM.normalize_by_col(psi_in_l, NR)
psi_in_r = psi_in_l
KPM.kpm_3d!(H_norm, Jx, Jy, Jz, NC, NR, NH, mu_3d_xyz, psi_in_l, psi_in_r; verbose=1, arr_size=16);

#fn = pt*"GammaxyzL100gam08D6NC512W01NR1Sam"*"$disind"*".jld2"
fn = "GammaxyzL100gam08D6NC512W2035NR1Sam"*"$disind"*".jld2"
jldsave(fn; Ls=Ls, NC=NC, D = D, Gammaxyz=mu_3d_xyz)

#=
mu_3d_yxz = zeros(ComplexF64, NC, NC, NC)
psi_in_l = exp.(2pi * 1im * rand(H_norm.n, NR));
KPM.normalize_by_col(psi_in_l, NR)
psi_in_r = psi_in_l
KPM.kpm_3d!(H_norm, Jy, Jx, Jz, NC, NR, NH, mu_3d_yxz, psi_in_l, psi_in_r; verbose=1, arr_size=16);

#fn = pt*"GammayxzL100gam08D6NC512W01NR1Sam"*"$disind"*".jld2"
fn = "GammayxzL100gam08D6NC512W2035NR1Sam"*"$disind"*".jld2"
jldsave(fn; Ls=Ls, NC=NC, D = D, Gammaxyz=mu_3d_yxz)
=#


