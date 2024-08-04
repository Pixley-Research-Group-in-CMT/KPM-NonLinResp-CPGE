using Statistics
using LinearAlgebra, Arpack, SparseArrays
using KPM
using JLD2
using CUDA

@show CUDA.has_cuda()
@show KPM.whichcore()
# check number of threads
println("Number of threads: $(Threads.nthreads())")

@show disind = parse(Int64,ARGS[1])
pt = "./CPGEdata/"  # path to save data
#pt = "/scratch/aw666/CPGEdata/"  # amarel
#fn = pt*"HdisL100gam08W2025dis"*"$disind"*".jld2"
#fn = pt*"HPertdisL100gam08W05PotDis"*"$disind"*".jld2"
fn = pt*"HdisL18gamma03W05PotDisTBC"*"$disind"*".jld2"
H_norm = load(fn, "H_norm")
Jx = load(fn, "Jx")
Jy = load(fn, "Jy")
Jz = load(fn, "Jz")

NC = 512
NR = 100
NH = H_norm.n
D = 6.0 #(M+2)+sqrt(1+gamma^2)
L = 18#100
Ls = [L,L,L]
#=
#fn = pt*"GammaxyzL100gam08D6NC512W2025NR1Sam"*"$disind"*".jld2"
fn = pt*"GammaxyzL$(L)gam03D6NC512W05NR$(NR)Sam"*"$disind"*".jld2"
@show isfile(fn)
if ! isfile(fn)
mu_3d_xyz = zeros(ComplexF64, NC, NC, NC)
psi_in_l = exp.(2pi * 1im * rand(H_norm.n, NR));
KPM.normalize_by_col(psi_in_l, NR)
psi_in_r = psi_in_l
KPM.kpm_3d!(H_norm, Jx, Jy, Jz, NC, NR, NH, mu_3d_xyz, psi_in_l, psi_in_r; verbose=1, arr_size=16);

#fn = pt*"GammaxyzL100gam08D6NC512W2005NR1Sam"*"$disind"*".jld2"
jldsave(fn; Ls=Ls, NC=NC, D = D, Gammaxyz=mu_3d_xyz)
end
=#

#fn = pt*"GammayxzL100gam08D6NC512W2025NR1Sam"*"$disind"*".jld2"
fn = pt*"GammayxzL$(L)gam03D6NC512W05NR$(NR)Sam"*"$disind"*".jld2"
if ! isfile(fn)
mu_3d_yxz = zeros(ComplexF64, NC, NC, NC)
psi_in_l = exp.(2pi * 1im * rand(H_norm.n, NR));
KPM.normalize_by_col(psi_in_l, NR)
psi_in_r = psi_in_l
KPM.kpm_3d!(H_norm, Jy, Jx, Jz, NC, NR, NH, mu_3d_yxz, psi_in_l, psi_in_r; verbose=1, arr_size=16);

#fn = pt*"GammayxzL100gam08D6NC512W2005NR1Sam"*"$disind"*".jld2"
jldsave(fn; Ls=Ls, NC=NC, D = D, Gammayxz=mu_3d_yxz)
end

